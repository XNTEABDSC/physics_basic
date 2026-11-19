use approx::AbsDiffEq;
use nalgebra::{Const, RealField, SVector};
use wacky_bag::utils::fn_as_index::FnAsIndex;
use std::{array, collections::HashSet, ops::{ControlFlow::Continue, Index, IndexMut}, path::Iter};

use crate::{rotation::DimSquare, shape::ShapeStruct};

/// 泛型凸多面体：由顶点和一组 (DIM-1) 维面构成。
/// 每个面是 DIM 个顶点确定的单纯形，`normal` 为外向单位法线，
/// `offset` 满足 `normal · x = offset` 对该面上所有点成立。
#[derive(Debug, Clone)]
pub struct ConvexPolyhedron<Num: RealField + Copy, const DIM: usize> {
    pub vertexes: Vec<SVector<Num, DIM>>,
    pub surfaces: Vec<Face<Num, DIM>>,
}

#[derive(Debug, Clone)]
pub struct Face<Num: RealField + Copy, const DIM: usize> {
    /// 该面的 DIM 个顶点索引。
    pub indices: [usize; DIM],
    /// 外向单位法线。
    pub normal: SVector<Num, DIM>,
    /// 平面偏移：normal · x = offset。
    pub offset: Num,
    /// 该面当前关联的外部点集合（仅 Quickhull 内部使用）。
    pub outside: Vec<usize>,
}

impl<Num: RealField + Copy, const DIM: usize> Face<Num, DIM> {
    fn signed_distance(&self, point: &SVector<Num, DIM>) -> Num {
        self.normal.dot(point) - self.offset
    }
}

impl<Num: RealField + Copy, const DIM: usize> ConvexPolyhedron<Num, DIM> {
    /// 从任意点集计算凸包。点数不足以构成 DIM 维单纯形时返回 `None`。
    pub fn quickhull(points: impl IntoIterator<Item = SVector<Num, DIM>>) -> Option<Self> {
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
            let (normal,offset) = make_outward_face(&mut indices,&self.vertexes, &centroid);
            self.surfaces.push(Face {
                indices,
                normal,
                offset,
                outside: Vec::new(),
            });
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
            // let verts: [SVector<Num, DIM>; DIM] = indices.map(|i| self.points[i]);
            let (normal,offset) = make_outward_face(&mut indices, &self.vertexes, &interior);
            // let offset = normal.dot(&verts[0]);
            self.surfaces.push(Face {
                indices,
                normal,
                offset,
                outside: Vec::new(),
            });
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

    // 构建 DIM-1 个差分向量
    let mut diffs = Vec::with_capacity(DIM - 1);
    for i in 1..DIM {
        diffs.push(verts[i] - verts[0]);
    }

    // 在 DIM 维中，DIM-1 个向量的法线可通过广义叉积得到：
    // normal_j = (-1)^j * det(矩阵去掉第 j 列)
    let mut normal: nalgebra::Matrix<Num, nalgebra::Const<DIM>, nalgebra::Const<1>, nalgebra::ArrayStorage<Num, DIM, 1>> = SVector::zeros();
    for j in 0..DIM {
        let mut sub = nalgebra::DMatrix::<Num>::zeros(DIM - 1, DIM - 1);
        let mut col = 0;
        for jj in 0..DIM {
            if jj == j {
                continue;
            }
            for i in 0..DIM - 1 {
                sub[(i, col)] = diffs[i][jj];
            }
            col += 1;
        }
        let det = sub.determinant();
        normal[j] = if j % 2 == 0 { det } else { -det };
    }

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
) -> (SVector<Num, DIM>, Num) {
    let verts: [SVector<Num, DIM>; DIM] = indices.map(|i| points[i]);
    let mut normal = face_normal(&verts);

    // 朝内 ⟺ 内部点落在法线的正侧
    if normal.dot(&(interior - verts[0])) > Num::zero() {
        indices.swap(0, 1);
        normal = -normal;
    }

    let offset = normal.dot(&verts[0]);
    (normal, offset)
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
        let hull = ConvexPolyhedron::quickhull(points).unwrap();

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
        let hull = ConvexPolyhedron::quickhull(points).unwrap();

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
            assert!(approx_eq(face.normal.norm(), 1.0));
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
        let hull = ConvexPolyhedron::quickhull(points).unwrap();

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
            assert!(approx_eq(face.normal.norm(), 1.0));
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
        assert!(ConvexPolyhedron::<f64, 2>::quickhull(points).is_none());
    }

    #[test]
    fn test_degenerate_coplanar_3d() {
        let points: Vec<SVector<f64, 3>> = vec![
            SVector::<f64, 3>::new(0.0, 0.0, 0.0),
            SVector::<f64, 3>::new(1.0, 0.0, 0.0),
            SVector::<f64, 3>::new(0.0, 1.0, 0.0),
            SVector::<f64, 3>::new(1.0, 1.0, 0.0),
        ];
        assert!(ConvexPolyhedron::<f64, 3>::quickhull(points).is_none());
    }
}