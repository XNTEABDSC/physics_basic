use approx::AbsDiffEq;
use nalgebra::{Const, DefaultAllocator, OMatrix, OVector, RealField, SVector};
use wacky_bag::utils::fn_as_index::FnAsIndex;
use std::{array, collections::HashSet, ops::{ControlFlow::Continue, Deref, Index, IndexMut}, path::Iter};

use crate::{shape::{self, ShapeStruct, generalized_cross, surface_center, volume}};
use crate::dimension::*;
/// 泛型凸多面体：由顶点和一组 (DIM-1) 维面构成。
/// 每个面是 DIM 个顶点确定的单纯形，`normal` 为外向单位法线，
/// `offset` 满足 `normal · x = offset` 对该面上所有点成立。
#[derive(Debug, Clone)]
pub struct ConvexPolyhedron<Num: RealField + Copy, const DIM: usize> {
    vertexes: Vec<SVector<Num, DIM>>,
    surfaces: Vec<Surface<Num, DIM>>,

    volume:Num
}

#[derive(Debug, Clone)]
pub struct Surface<Num: RealField + Copy, const DIM: usize> {
    /// 该面的 DIM 个顶点索引。
    indices: [usize; DIM],
    /// 外向单位法线。
    normal_vector: SVector<Num, DIM>,
    /// 平面偏移：normal · x = offset。
    offset: Num,
    /// 该面当前关联的外部点集合（仅 Quickhull 内部使用）。
    outside: Vec<usize>,

    /// 面向量，方向为法向量，长度为面积
	edge_vector:SVector<Num,DIM>,
    area:Num,
	center:SVector<Num,DIM>
}

impl<Num: RealField + Copy, const DIM: usize> Surface<Num, DIM> {
    fn signed_distance(&self, point: &SVector<Num, DIM>) -> Num {
        self.normal_vector.dot(point) - self.offset
    }

    fn new(
        mut indices:[usize; DIM],
        vertexes: &[SVector<Num, DIM>],
        interior: &SVector<Num, DIM>,
    )->Self{
        let (edge_vector, normal_vector,offset) = make_outward_face(&mut indices,vertexes, interior);
        Surface {
            indices,
            normal_vector: normal_vector,
            offset,
            outside: Vec::new(),

            area: edge_vector.magnitude(),
            edge_vector: edge_vector,

            center:surface_center(vertexes,&indices)
        }
    }
}

impl<Num: RealField + Copy, const DIM: usize> ConvexPolyhedron<Num, DIM> 
where Const<DIM>:DimSquare
{
    /// 从任意点集计算凸包。点数不足以构成 DIM 维单纯形时返回 `None`。
    pub fn new(points: impl IntoIterator<Item = SVector<Num, DIM>>) -> Option<Self> {
        assert!(DIM >= 2, "DIM 必须 >= 2");

        // 1. 去重
        let mut unique = Vec::<SVector<Num, DIM>>::new();
        // let mut seen = HashSet::new();
        for p in points {
            if unique.iter().all(|v| !v.abs_diff_eq(&p, Num::default_epsilon()))
            {
                unique.push(p);
            }
        }
        if unique.len() < DIM + 1 {
            return None;
        }

        // 2. 寻找初始单纯形（DIM+1 个仿射无关的点）
        let simplex = find_initial_simplex(&unique)?;
        let mut hull = Self {
            vertexes: unique,
            surfaces: Vec::new(),

            volume:Num::zero()
        };
        hull.build_initial_faces(&simplex);

        // 3. 将所有剩余点分配到可见面
        let mut unassigned: Vec<usize> = (0..hull.vertexes.len())
            .filter(|i| !simplex.contains(i))
            .collect();

        for &pi in &unassigned {
            let p = hull.vertexes[pi];
            for face in &mut hull.surfaces {
                if face.signed_distance(&p) > Num::default_epsilon() {
                    face.outside.push(pi);
                }
            }
        }

        // 4. Quickhull 主循环
        loop {
            // 找到拥有最远外部点的面
            let mut best_face = None;
            let mut best_dist = Num::zero();
            for (fi, face) in hull.surfaces.iter().enumerate() {
                if let Some(&pi) = face
                    .outside
                    .iter()
                    .max_by(|&&a, &&b| {
                        face.signed_distance(&hull.vertexes[a])
                            .partial_cmp(&face.signed_distance(&hull.vertexes[b]))
                            .unwrap()
                    })
                {
                    let d = face.signed_distance(&hull.vertexes[pi]);
                    if d > best_dist {
                        best_dist = d;
                        best_face = Some((fi, pi));
                    }
                }
            }

            let Some((start_face, apex)) = best_face else {
                break;
            };

            // 删除所有可见面，收集地平线
            let horizon = hull.remove_visible_faces(apex);

            // 用地平线上的每条边与 apex 构成新面
            hull.create_new_faces(apex, &horizon);

            // 将被删除面的外部点重新分配给新面
            // （简化版：重新扫描所有未归属点；实际实现可用可见面的外部点集合提高效率）
            let assigned: HashSet<usize> = hull
                .surfaces
                .iter()
                .flat_map(|f| f.outside.iter().copied())
                .collect();
            for pi in 0..hull.vertexes.len() {
                if pi == apex || simplex.contains(&pi) || assigned.contains(&pi) {
                    continue;
                }
                let p = hull.vertexes[pi];
                for face in &mut hull.surfaces {
                    if face.signed_distance(&p) > Num::default_epsilon()
                        && !face.outside.contains(&pi)
                    {
                        face.outside.push(pi);
                    }
                }
            }
        }

        // 清理：移除没有顶点的面（理论上不会出现）
        hull.surfaces.retain(|f| f.indices.len() == DIM);

        hull.volume=volume(&hull.vertexes, hull.surfaces.iter().map(|v|&v.indices));
        Some(hull)
    }

    fn build_initial_faces(&mut self, simplex: &ArrayChain<usize,DIM,1>) {
        // 对每个 DIM 元子集构建一个面，法线方向由单纯形的重心决定
        let centroid = simplex
            .iter()
            .map(|&i| self.vertexes[i])
            .fold(SVector::zeros(), |acc, p| acc + p)
            / Num::from_usize(simplex.len()).unwrap();

        for k in 0..=DIM {
            let mut idx = Vec::with_capacity(DIM);
            for j in 0..=DIM {
                if j != k {
                    idx.push(simplex[j]);
                }
            }
            let mut indices = [0usize; DIM];
            indices.copy_from_slice(&idx);

            let surface=Surface::new(indices, &self.vertexes, &centroid);
            self.surfaces.push(surface);
            // let (normal,offset) = make_outward_face(&mut indices,&self.vertexes, &centroid);
            // let normal=normal.
            // self.surfaces.push(Surface {
            //     indices,
            //     normal_vector: normal,
            //     offset,
            //     outside: Vec::new(),
            // });
        }
    }

    fn remove_visible_faces(&mut self, apex: usize) -> HashSet<Vec<usize>> {
        let apex_pos = self.vertexes[apex];
        let mut visible = Vec::new();
        for (i, face) in self.surfaces.iter().enumerate() {
            if face.signed_distance(&apex_pos) > Num::default_epsilon() {
                visible.push(i);
            }
        }

        let mut horizon = HashSet::new();
        for &fi in &visible {
            let face = &self.surfaces[fi];
            // 一个 face 是 (DIM-1)-单纯形（DIM 个顶点），
            // 其边界 facet 是去掉某个顶点得到的 (DIM-2)-单纯形（DIM-1 个顶点）
            for k in 0..DIM {
                let mut sub: Vec<usize> =
                    (0..DIM).filter(|&j| j != k).map(|j| face.indices[j]).collect();
                sub.sort_unstable();
                if !horizon.insert(sub.clone()) {
                    horizon.remove(&sub);      // 两个可见面共享 → 内部，不是地平线
                }
            }
        }

        visible.sort_unstable_by(|&a, &b| b.cmp(&a));
        for &fi in &visible {
            self.surfaces.swap_remove(fi);
        }
        horizon
    }

    fn create_new_faces(&mut self, apex: usize, horizon: &HashSet<Vec<usize>>) {
        let interior = self.interior_point();
        for sub in horizon {
            debug_assert_eq!(sub.len(), DIM - 1);
            let mut indices = [0usize; DIM];
            indices[0] = apex;
            for (i, &v) in sub.iter().enumerate() {
                indices[i + 1] = v;
            }
            
            let surface=Surface::new(indices, &self.vertexes, &interior);
            self.surfaces.push(surface);
            // let verts: [SVector<Num, DIM>; DIM] = indices.map(|i| self.points[i]);
            // let (normal,offset) = make_outward_face(&mut indices, &self.vertexes, &interior);
            // // let offset = normal.dot(&verts[0]);
            // self.surfaces.push(Surface {
            //     indices,
            //     normal_vector: normal,
            //     offset,
            //     outside: Vec::new(),
            // });
        }
    }

    fn interior_point(&self) -> SVector<Num, DIM> {
        let n = self.vertexes.len();
        self.vertexes.iter().fold(SVector::zeros(), |acc, p| acc + p)
            / Num::from_usize(n).unwrap()
    }

    pub fn into_shape_poly(self)->ShapeStruct<Num,DIM>
    where Const<DIM>:DimSquare
    {
        ShapeStruct::new(self.vertexes, self.surfaces.into_iter().map(|i|i.indices).collect())
    }


    /// 数值容差，取机器精度的平方根，与标量类型精度相关。
    fn eps() -> Num {
        Num::default_epsilon().sqrt()
    }

    /// 判断点是否在凸多面体内部（含边界）。
    pub fn contains_point(&self, p: &SVector<Num, DIM>) -> bool {
        let e = Self::eps();
        for s in &self.surfaces {
            // 外向法线：内部满足 normal · x <= offset
            if s.normal_vector.dot(p) - s.offset > e {
                return false;
            }
        }
        true
    }

    /// 计算直线 origin + t * direction 与多面体的交集。
    /// 返回相交线段的两个端点；若无交返回 None。
    /// 若直线与多面体相切于一点，返回两个相同的点。
    /// ```
    /// let p0 = origin + direction * tmin;
    /// let p1 = origin + direction * tmax;
    /// ```
    pub fn line_intersection(
        &self,
        origin: &SVector<Num, DIM>,
        direction: &SVector<Num, DIM>,
    ) -> Option<(Num,Num)> {
        let e = Self::eps();

        // 退化为点
        if direction.norm_squared() < e * e {
            return if self.contains_point(origin) {
                Some((Num::zero(),Num::zero()))
            } else {
                None
            };
        }

        let mut tmin: Option<Num> = None;
        let mut tmax: Option<Num> = None;

        for s in &self.surfaces {
            // 约束：normal · (origin + t * direction) <= offset
            // 即：t * (normal · direction) <= offset - normal · origin
            let a = s.normal_vector.dot(direction);
            let b = s.offset - s.normal_vector.dot(origin);

            if a.abs() < e {
                // 直线方向与该面平行
                if b < -e {
                    // 整条直线都在该面外侧
                    return None;
                }
            } else {
                let t = b / a;
                if a > Num::zero() {
                    // t <= t
                    if tmax.map_or(true, |tm| t < tm) {
                        tmax = Some(t);
                    }
                } else {
                    // t >= t
                    if tmin.map_or(true, |tm| t > tm) {
                        tmin = Some(t);
                    }
                }
            }
        }

        // 凸多面体有界，tmin / tmax 应均存在
        let tmin = tmin?;
        let tmax = tmax?;

        if tmin > tmax {
            // 允许相切数值误差
            if tmin - tmax < e {
                let t = (tmin + tmax) / Num::from_usize(2).unwrap();
                // let p = origin + direction * t;
                return Some((t, t));
            }
            return None;
        }

        // let p0 = origin + direction * tmin;
        // let p1 = origin + direction * tmax;
        Some((tmin, tmax))
    }

    /// 顶点质心（作为 GJK 初始方向）。
    fn centroid(&self) -> SVector<Num, DIM> {
        let mut c = SVector::<Num, DIM>::zeros();
        if self.vertexes.is_empty() {
            return c;
        }
        for v in &self.vertexes {
            c += v;
        }
        c / Num::from_usize(self.vertexes.len()).unwrap()
    }

    /// 在给定方向上求支撑点（点积最大的顶点）。
    fn support(&self, dir: &SVector<Num, DIM>) -> SVector<Num, DIM> {
        let mut best = self.vertexes[0];
        let mut best_dot = dir.dot(&best);
        for v in &self.vertexes[1..] {
            let d = dir.dot(v);
            if d > best_dot {
                best_dot = d;
                best = *v;
            }
        }
        best
    }

    /// 两个凸多面体是否相交（碰撞检测）。
    /// 使用 GJK 判断 Minkowski 差 A - B 是否包含原点。
    pub fn intersects(&self, other: &Self) -> bool 
    where 
        Const<DIM>:DimToDimAdd1,
        DefaultAllocator:AllocatorVM<ConstDimToDimAdd1T<DIM>>,
        ConstDimToDimAdd1T<DIM>:DimSquare
    {
        let e = Self::eps();

        // 初始搜索方向
        let mut dir = self.centroid() - other.centroid();
        if dir.norm_squared() < e * e {
            dir = SVector::<Num, DIM>::zeros();
            dir[0] = Num::one();
        }

        // Minkowski 差 A - B 的支撑函数
        let support_diff = |d: &SVector<Num, DIM>| -> SVector<Num, DIM> {
            self.support(d) - other.support(&(-d))
        };

        let first = support_diff(&dir);
        let mut simplex: Vec<SVector<Num, DIM>> = vec![first];
        let mut v = -first; // 指向原点的方向

        if v.norm_squared() < e * e {
            return true; // 原点已在 Minkowski 差中
        }

        loop {
            let w = support_diff(&v);

            // w 不能使单纯形在 -v 方向上越过原点，则不相交
            if v.dot(&w) < -e {
                return false;
            }

            // 支撑点重复，无法继续推进
            if simplex.iter().any(|p| (*p - w).norm_squared() < e * e) {
                let (cp, _) = Self::closest_point_on_simplex(&simplex);
                return cp.norm_squared() < e * e;
            }

            simplex.push(w);
            let (cp, subset) = Self::closest_point_on_simplex(&simplex);

            if cp.norm_squared() < e * e {
                return true; // 原点在单纯形内 ⇒ Minkowski 差包含原点
            }

            // 缩减单纯形到包含最近点的最小面
            let mut new_simplex = Vec::with_capacity(subset.len());
            for &i in &subset {
                new_simplex.push(simplex[i]);
            }
            simplex = new_simplex;
            v = -cp;
        }
    }

    /// 求一组点（单纯形）的凸包上离原点最近的点，并返回其支撑子集。
    /// 枚举所有非空子集，对每个子集解最小二乘问题：
    ///   min ||Σ λ_i p_i||²  s.t.  Σ λ_i = 1
    /// 取满足 λ_i ≥ 0 且范数最小的解。
    fn closest_point_on_simplex(
        points: &[SVector<Num, DIM>],
    ) -> (SVector<Num, DIM>, Vec<usize>) 
    where 
        Const<DIM>:DimToDimAdd1,
        DefaultAllocator:AllocatorVM<ConstDimToDimAdd1T<DIM>>,
        ConstDimToDimAdd1T<DIM>:DimSquare
    {
        let e = Self::eps();
        let n = points.len();
        let mut best_norm2: Option<Num> = None;
        let mut best_point = SVector::<Num, DIM>::zeros();
        let mut best_subset: Vec<usize> = vec![];

        // 枚举所有非空子集
        for mask in 1..(1usize << n) {
            let subset: Vec<usize> = (0..n).filter(|i| mask & (1 << i) != 0).collect();
            let k = subset.len();

            // 构造 KKT 系统：
            // [ G  1 ] [λ]   [0]
            // [ 1ᵀ 0 ] [μ] = [1]
            // 其中 G_{ij} = p_i · p_j
            let mut mat = OMatrix::<Num,ConstDimToDimAdd1T<DIM>,ConstDimToDimAdd1T<DIM>>::zeros();
            for i in 0..k {
                for j in 0..k {
                    mat[(i, j)] = points[subset[i]].dot(&points[subset[j]]);
                }
                mat[(i, k)] = Num::one();
                mat[(k, i)] = Num::one();
            }
            let mut rhs = OVector::<Num,ConstDimToDimAdd1T<DIM>>::zeros();
            rhs[k] = Num::one();

            if let Some(sol) = mat.lu().solve(&rhs) {
                // 检查重心坐标非负
                let mut ok = true;
                for i in 0..k {
                    if sol[i] < -e {
                        ok = false;
                        break;
                    }
                }
                if ok {
                    let mut pt = SVector::<Num, DIM>::zeros();
                    for i in 0..k {
                        pt += points[subset[i]] * sol[i];
                    }
                    let norm2 = pt.norm_squared();
                    if best_norm2.map_or(true, |b| norm2 < b) {
                        best_norm2 = Some(norm2);
                        best_point = pt;
                        best_subset = subset;
                    }
                }
            }
        }

        (best_point, best_subset)
    }
}

impl<Num: RealField + Copy, const DIM: usize> shape::Surface<Num,DIM> for  Surface<Num, DIM>  {
    fn vertexes_index(&self)->impl std::ops::Deref<Target = [usize;DIM]> {
        &self.indices
    }

    fn edge_vector(&self)->impl std::ops::Deref<Target = SVector<Num,DIM>> {
        &self.edge_vector
    }

    fn normal_vector(&self)->impl std::ops::Deref<Target = SVector<Num,DIM>> {
        &self.normal_vector
    }

    fn area(&self)->impl std::ops::Deref<Target = Num> {
        &self.area
    }

    fn center(&self)->impl std::ops::Deref<Target = SVector<Num,DIM>> {
        &self.center
    }
}

impl<Num: RealField + Copy, const DIM: usize> shape::Shape<Num,DIM> for ConvexPolyhedron<Num,DIM>  {
    fn surfaces(&self)->impl
    Deref<Target =Self::SurfaceCollection> {
        &self.surfaces
    }

    fn vertexes(&self)->impl std::ops::Deref<Target : Index<usize,Output = SVector<Num,DIM>>> {
        &self.vertexes
    }

    fn volume(&self)->impl std::ops::Deref<Target = Num> {
        &self.volume
    }
    
    type Surface=Surface<Num,DIM>;
    
    type SurfaceCollection=Vec<Surface<Num, DIM>>;
}

struct ArrayChain<T,const N1:usize,const N2:usize>(pub [T;N1],pub [T;N2]);

impl<T,const N1:usize,const N2:usize> ArrayChain<T,N1,N2>{
    pub fn len(&self)->usize{N1+N2}
    pub fn from_fn<F>(mut f:F)->Self
    where F:FnMut(usize)->T
    {
        Self(
            array::from_fn(&mut f),
            array::from_fn(move |i|f(i+N1)) 
        )
    }
    pub fn contains(&self,v:&T)->bool
    where
        T: PartialEq,
    {
        return self.0.contains(v) || self.1.contains(v);
    }
    pub fn iter(&self)->impl Iterator<Item=&T>{
        self.0.iter().chain(self.1.iter())
    }
}

impl<T,const N1:usize,const N2:usize> Index<usize> for ArrayChain<T,N1,N2> {
    type Output=T;

    fn index(&self, index: usize) -> &Self::Output {
        if index<N1 {
            &self.0[index]
        } else {
            &self.1[index-N1]
        }
    }
}
impl<T,const N1:usize,const N2:usize> IndexMut<usize> for ArrayChain<T,N1,N2> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        if index<N1 {
            &mut self.0[index]
        } else {
            &mut self.1[index-N1]
        }
    }
}


/// 寻找 DIM+1 个仿射无关的点。
fn find_initial_simplex<Num: RealField + Copy, const DIM: usize>(
    points: &[SVector<Num, DIM>],
) -> Option<ArrayChain<usize,DIM,1>> {
    let n = points.len();

    let mut chosen=ArrayChain::from_fn(|_|n+1);
    // // 第一个点：任意
    chosen[0] = 0;

    for k in 1..=DIM {
        let base = &chosen;
        let mut found = None;
        for i in 0..n {
            if base.contains(&i) {
                continue;
            }
            // 检查是否与已有 k 个点仿射无关
            if affine_independent( FnAsIndex::new(|i|&points[base[i]]), &points[i], k) {
                found = Some(i);
                break;
            }
        }
        chosen[k]=found?;
    }
    Some(chosen)
}

fn affine_independent<Num: RealField + Copy, const DIM: usize>(
    base: impl Index<usize,Output = SVector<Num, DIM>>,
    candidate: &SVector<Num, DIM>,
    k: usize,
) -> bool {
    if k == 0 {
        return true;
    }
    // 在 k 维仿射空间中，候选点不能落在 base 张成的仿射子空间中。
    // 简化检查：到 base 的仿射包的距离 > eps。
    // 对于 k=1: 不能等于 base[0]
    // k=2: 不能共线
    // k=3: 不能共面
    // 通用做法：构建差分矩阵，检查秩。
    let mut diffs = Vec::with_capacity(k);
    for i in 1..k {
        diffs.push(base[i] - base[0]);
    }
    diffs.push(*candidate - base[0]);
    let rank = rank_of_vectors(&diffs);
    rank == k
}

fn rank_of_vectors<Num: RealField + Copy, const DIM: usize>(
    vectors: &[SVector<Num, DIM>],
) -> usize {
    let mut mat = nalgebra::DMatrix::<Num>::zeros(vectors.len(), DIM);
    for (i, v) in vectors.iter().enumerate() {
        for j in 0..DIM {
            mat[(i, j)] = v[j];
        }
    }
    mat.rank(Num::default_epsilon())
}

/// 计算 DIM 个顶点确定的超平面的外向单位法线。
/// `interior` 是凸包内部一点，用于确定法线方向。
fn face_normal<Num: RealField + Copy, const DIM: usize>(
    verts: &[SVector<Num, DIM>; DIM],
) -> SVector<Num, DIM> {
    assert_eq!(verts.len(), DIM);

    let mut normal=generalized_cross(verts.clone());
    // // 构建 DIM-1 个差分向量
    // let mut diffs = Vec::with_capacity(DIM - 1);
    // for i in 1..DIM {
    //     diffs.push(verts[i] - verts[0]);
    // }

    // // 在 DIM 维中，DIM-1 个向量的法线可通过广义叉积得到：
    // // normal_j = (-1)^j * det(矩阵去掉第 j 列)
    // let mut normal: nalgebra::Matrix<Num, nalgebra::Const<DIM>, nalgebra::Const<1>, nalgebra::ArrayStorage<Num, DIM, 1>> = SVector::zeros();
    // for j in 0..DIM {
    //     let mut sub = nalgebra::DMatrix::<Num>::zeros(DIM - 1, DIM - 1);
    //     let mut col = 0;
    //     for jj in 0..DIM {
    //         if jj == j {
    //             continue;
    //         }
    //         for i in 0..DIM - 1 {
    //             sub[(i, col)] = diffs[i][jj];
    //         }
    //         col += 1;
    //     }
    //     let det = sub.determinant();
    //     normal[j] = if j % 2 == 0 { det } else { -det };
    // }

    if normal.try_normalize_mut(Num::default_epsilon()).is_none(){normal=SVector::zeros();}

    normal
}

/// 保证 (indices 顺序 → 法线) 这条约定成立，且法线朝外。
/// 若按 indices 顺序得到的法线指向内部，则 `indices.swap(0, 1)` 翻转顺序，
/// 从而翻转法线。返回最终的 (外向单位法线, offset)。
fn make_outward_face<Num: RealField + Copy, const DIM: usize>(
    indices: &mut [usize; DIM],
    points: &[SVector<Num, DIM>],
    interior: &SVector<Num, DIM>,
) -> (SVector<Num, DIM>, SVector<Num, DIM>, Num) {
    let verts: [SVector<Num, DIM>; DIM] = indices.map(|i| points[i]);
    let mut edge_vector = generalized_cross(verts);

    // 朝内 ⟺ 内部点落在法线的正侧
    if edge_vector.dot(&(interior - verts[0])) > Num::zero() {
        indices.swap(0, 1);
        edge_vector = -edge_vector;
    }
    let mut normal=edge_vector;
    if normal.try_normalize_mut(Num::default_epsilon()).is_none(){normal=SVector::zeros();}
    let offset = normal.dot(&verts[0]);
    (edge_vector, normal, offset)
}

/// 在 DIM>3 时，用离 apex、a、b 仿射无关的点补齐面索引。
fn fill_face_indices<Num: RealField + Copy, const DIM: usize>(
    points: &[SVector<Num, DIM>],
    indices: &mut [usize; DIM],
    apex: usize,
    a: usize,
    b: usize,
) {
    let mut filled = 3;
    for i in 0..points.len() {
        if filled == DIM {
            break;
        }
        if i == apex || i == a || i == b {
            continue;
        }
        let base: Vec<_> = indices[0..filled].iter().map(|&j| points[j]).collect();
        if affine_independent(base, &points[i], filled) {
            indices[filled] = i;
            filled += 1;
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::SVector;

    fn approx_eq(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn test_2d_square() {
        // 正方形 + 内部点
        let points: Vec<SVector<f64, 2>> = vec![
            SVector::<f64, 2>::new(0.0, 0.0),
            SVector::<f64, 2>::new(1.0, 0.0),
            SVector::<f64, 2>::new(1.0, 1.0),
            SVector::<f64, 2>::new(0.0, 1.0),
            SVector::<f64, 2>::new(0.5, 0.5),
        ];
        let hull = ConvexPolyhedron::new(points).unwrap();

        // println!("{:?}",hull);

        // 2D 凸包的面是边，每个面由 2 个顶点构成
        for face in &hull.surfaces {
            assert_eq!(face.indices.len(), 2);
        }

        // 顶点数：4 个角点 + 可能包含内部点（取决于实现是否保留内部点）
        // 我们的实现保留所有去重后的点，但内部点不会出现在任何面的 indices 中
        let hull_vertex_indices: HashSet<usize> =
            hull.surfaces.iter().flat_map(|f| f.indices).collect();

        // 每个顶点都在某条边上
        for i in 0..4 {
            assert!(hull_vertex_indices.contains(&i),"{i} should in hull"); // error 3 should in hull
        }
        // 内部点 (0.5, 0.5) 不应出现在任何面上
        assert!(!hull_vertex_indices.contains(&4),"{} should not in hull",4);
    }

    #[test]
    fn test_3d_cube() {
        // 立方体 8 个顶点 + 中心点
        let points: Vec<SVector<f64, 3>> = vec![
            SVector::<f64, 3>::new(0.0, 0.0, 0.0),
            SVector::<f64, 3>::new(1.0, 0.0, 0.0),
            SVector::<f64, 3>::new(1.0, 1.0, 0.0),
            SVector::<f64, 3>::new(0.0, 1.0, 0.0),
            SVector::<f64, 3>::new(0.0, 0.0, 1.0),
            SVector::<f64, 3>::new(1.0, 0.0, 1.0),
            SVector::<f64, 3>::new(1.0, 1.0, 1.0),
            SVector::<f64, 3>::new(0.0, 1.0, 1.0),
            SVector::<f64, 3>::new(0.5, 0.5, 0.5),
        ];
        let hull = ConvexPolyhedron::new(points).unwrap();

        // 3D 凸包的面是三角形
        for face in &hull.surfaces {
            assert_eq!(face.indices.len(), 3);
        }

        // 立方体凸包应有 8 个顶点参与
        let hull_vertices: HashSet<usize> =
            hull.surfaces.iter().flat_map(|f| f.indices).collect();
        assert_eq!(hull_vertices.len(), 8);

        // 每个面的法线应该是单位向量，且所有立方体顶点都在对应面的内部或负侧
        for face in &hull.surfaces {
            assert!(approx_eq(face.normal_vector.norm(), 1.0));
            for (i, &pi) in hull_vertices.iter().enumerate() {
                let d = face.signed_distance(&hull.vertexes[pi]);
                // 凸包上的点应满足 signed_distance <= eps
                assert!(d <= 1e-9, "face {} violates vertex {}", face.indices[0], i);
            }
        }
    }

    #[test]
    fn test_4d_hypercube() {
        // 4D 超立方体 16 个顶点
        let mut points = Vec::new();
        for i in 0..16 {
            let mut p = SVector::<f64, 4>::zeros();
            for d in 0..4 {
                p[d] = if (i >> d) & 1 == 0 { 0.0 } else { 1.0 };
            }
            points.push(p);
        }
        let hull = ConvexPolyhedron::new(points).unwrap();

        // 4D 凸包的面是四面体（4 个顶点）
        for face in &hull.surfaces {
            assert_eq!(face.indices.len(), 4);
        }

        // 超立方体凸包应包含全部 16 个顶点
        let hull_vertices: HashSet<usize> =
            hull.surfaces.iter().flat_map(|f| f.indices).collect();
        assert_eq!(hull_vertices.len(), 16);

        // 所有顶点都应满足所有面的半空间约束
        for face in &hull.surfaces {
            assert!(approx_eq(face.normal_vector.norm(), 1.0));
            for &pi in &hull_vertices {
                let d = face.signed_distance(&hull.vertexes[pi]);
                assert!(d <= 1e-9);
            }
        }
    }

    #[test]
    fn test_degenerate_collinear_2d() {
        let points: Vec<SVector<f64, 2>> = vec![
            SVector::<f64, 2>::new(0.0, 0.0),
            SVector::<f64, 2>::new(1.0, 1.0),
            SVector::<f64, 2>::new(2.0, 2.0),
        ];
        assert!(ConvexPolyhedron::<f64, 2>::new(points).is_none());
    }

    #[test]
    fn test_degenerate_coplanar_3d() {
        let points: Vec<SVector<f64, 3>> = vec![
            SVector::<f64, 3>::new(0.0, 0.0, 0.0),
            SVector::<f64, 3>::new(1.0, 0.0, 0.0),
            SVector::<f64, 3>::new(0.0, 1.0, 0.0),
            SVector::<f64, 3>::new(1.0, 1.0, 0.0),
        ];
        assert!(ConvexPolyhedron::<f64, 3>::new(points).is_none());
    }
}