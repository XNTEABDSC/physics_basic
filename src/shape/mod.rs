pub mod convex_polyhedron;
pub mod multi_convex;

use std::{array, ops::{Deref, Index}, process::Output};

use nalgebra::{Const, RealField, SMatrix, SVector};
use wacky_bag::{traits::as_iter::AsIterator, utils::{d_sphere_volume::d_sphere_volume_by_radius_pow, factorial::factorial}};
use crate::dimension::*;
use crate::{stats::Volume};

#[derive(Default, Clone, Copy, Debug)]
pub struct ShapeSphere<Num:RealField,const DIM:usize>{
    pub radius:Num,
    pub radius_pow:Num,
	pub volume:Num
}

impl<Num:RealField+Copy+num_traits::Num,const DIM:usize> ShapeSphere<Num,DIM> {
    pub fn from_radius(radius:Num)->Self {
		let radius_pow=radius.powi(DIM as i32);
        return Self{
            radius,
            radius_pow,
			volume:d_sphere_volume_by_radius_pow(radius_pow, DIM)
        };
    }

	pub fn volume(&self)->Volume<Num>{Volume(self.volume)}
}

#[derive(Debug,Clone)]
pub struct SurfaceStruct<Num:RealField,const DIM:usize>{
	pub vertexes_index:[usize;DIM],
	pub edge_vector:SVector<Num,DIM>,
    pub normal_vector:SVector<Num,DIM>,
    pub area:Num,
	pub center:SVector<Num,DIM>
}
#[derive(Debug,Clone)]
pub struct ShapeStruct<Num:RealField,const DIM:usize>{
	pub vertexes:Vec<SVector<Num,DIM>>,
	pub surfaces:Vec<SurfaceStruct<Num,DIM>>,
	pub volume:Num
}

pub trait Shape<Num,const DIM:usize> {
    // type Num;
    // const DIM:usize;

    type Surface:Surface<Num,DIM>;
    type SurfaceCollection:Index<usize,Output = Self::Surface>+for<'a> AsIterator< Item<'a> =&'a Self::Surface>;
    
    fn surfaces(&self)->impl Deref<Target = Self::SurfaceCollection>;
    fn vertexes(&self)->impl Deref<Target : Index<usize,Output = SVector<Num,DIM>>>;
    fn volume(&self)->impl Deref<Target = Num>;
}

impl<Num:RealField,const DIM:usize> Shape<Num,DIM> for ShapeStruct<Num,DIM>  {

    fn surfaces(&self)->impl Deref<Target=Self::SurfaceCollection>
    {
        &self.surfaces
    }

    fn vertexes(&self)->impl Deref<Target : Index<usize,Output = SVector<Num,DIM>>> {
        &self.vertexes
    }

    fn volume(&self)->impl Deref<Target = Num> {
        &self.volume
    }
    
    type Surface=SurfaceStruct<Num,DIM>;
    
    type SurfaceCollection=Vec<SurfaceStruct<Num, DIM>>;
    
}

pub trait Surface<Num,const DIM:usize> {
    
	fn vertexes_index(&self)->impl Deref<Target = [usize;DIM]>;
	fn edge_vector(&self)->impl Deref<Target = SVector<Num,DIM>>;
    fn normal_vector(&self)->impl Deref<Target = SVector<Num,DIM>>;
    fn area(&self)->impl Deref<Target = Num>;
	fn center(&self)->impl Deref<Target = SVector<Num,DIM>>;
}

impl<Num:RealField,const DIM:usize> Surface<Num,DIM> for SurfaceStruct<Num,DIM>  {
    fn vertexes_index(&self)->impl Deref<Target = [usize;DIM]> {
        &self.vertexes_index
    }

    fn edge_vector(&self)->impl Deref<Target = SVector<Num,DIM>> {
        &self.edge_vector
    }

    fn normal_vector(&self)->impl Deref<Target = SVector<Num,DIM>> {
        &self.normal_vector
    }

    fn area(&self)->impl Deref<Target = Num> {
        &self.area
    }

    fn center(&self)->impl Deref<Target = SVector<Num,DIM>> {
        &self.center
    }
}
/// 计算一个面的“未归一化”法向量。
/// DIM=3 时就是三角形叉积；
/// DIM=2 时就是线段的垂直向量；
/// 一般 DIM 下使用余子式/广义叉积。
pub fn generalized_cross<Num: RealField, const DIM: usize>(
    vertexes: impl IntoIterator<Item = SVector<Num,DIM>>,
    // edge: &[usize; DIM],
) -> SVector<Num, DIM> {
    assert!(DIM >= 2);
    let mut vertexes=vertexes.into_iter();
    let p0 = vertexes.next().unwrap();

    // DIM 维空间中，(DIM-1)-单纯形有 DIM 个顶点，
    // 因此有 DIM-1 条边向量。
    let m = DIM - 1;
    let mut edges = Vec::with_capacity(m);

    for i in 0..m {
        edges.push(vertexes.next().unwrap() - p0.clone());
    }

    // n[k] = (-1)^k * det(去掉第 k 列的 (DIM-1)x(DIM) 边矩阵)
    let mut n = SVector::<Num, DIM>::zeros();

    for k in 0..DIM {
        let mut mat = nalgebra::DMatrix::<Num>::zeros(m, m);

        for i in 0..m {
            let mut col = 0;
            for j in 0..DIM {
                if j == k {
                    continue;
                }
                mat[(i, col)] = edges[i][j].clone();
                col += 1;
            }
        }

        let mut d = mat.determinant();
        if k % 2 == 1 {
            d = -d;
        }
        n[k] = d;
    }

    n
}

pub fn surface_center<Num: RealField, const DIM: usize>(
    vertexes: &[SVector<Num, DIM>],
    edge: &[usize; DIM],
) -> SVector<Num, DIM> 
{
	let a:SMatrix<Num,DIM,DIM>=SMatrix::from( array::from_fn::<_,DIM,_>(|i|vertexes[edge[i]].clone().into())  );
	a.column_sum().scale(Num::one()/Num::from_isize(DIM as isize).unwrap())
}

impl<Num:RealField,const DIM:usize> SurfaceStruct<Num,DIM> {
	pub fn new(vertexes:&[SVector<Num,DIM>],edge:[usize;DIM])->Self{
		let edge_vector=generalized_cross(edge.iter().map(|i|vertexes[*i].clone()));
		let center=surface_center(vertexes,&edge);
		Self{
			vertexes_index:edge,
			normal_vector:edge_vector.normalize(),
			area:edge_vector.magnitude(),
			edge_vector,
			center
		}
	}
}

pub fn volume<'a,Num:RealField,const DIM:usize>(vertexes:&'a [SVector<Num,DIM>],edges:impl IntoIterator<Item = &'a [usize;DIM]>) -> Num 
where Const<DIM>:DimSquare
{
    let mut sum = Num::zero();

    for edges in edges {
        let mat = SMatrix::<Num, DIM, DIM>::from_fn(|r, c| {
            vertexes[edges[r]][c].clone()
        });

        sum += mat.determinant();
    }

    let fact = Num::from_usize(factorial(DIM).unwrap()).unwrap();

    (sum / fact).abs()
}

impl<Num:RealField,const DIM:usize> ShapeStruct<Num,DIM>
where Const<DIM>:DimSquare
{
	pub fn new(vertexes:Vec<SVector<Num,DIM>>,edges:Vec<[usize;DIM]>)->Self{
		let surfaces=edges.into_iter().map(|e|SurfaceStruct::new(&vertexes,e)).collect::<Vec<_>>();
		Self{
			volume:volume(&vertexes,surfaces.iter().map(|a|&a.vertexes_index)),
			vertexes,
			surfaces,
		}
	}
	

	/// 生成 N 维球体的多面体近似。
	///
	/// 参数：
	/// - `center`：球心
	/// - `radius`：半径
	/// - `segments`：每条轴上把 [-1,1] 分成的区间数（>=1）。
	///   `segments == 1` 即"超立方体角点投影"版本。
	///
	/// 注意：面数约为 `2·DIM·segments^(DIM-1)·(DIM-1)!`，
	/// 大 DIM 时建议只用 `segments = 1`。
	pub fn from_ball_grid(
		center: SVector<Num, DIM>,
		radius: Num,
		segments: usize,
	) -> ShapeStruct<Num, DIM> {
		assert!(DIM >= 2, "DIM must be >= 2");
		assert!(segments >= 1, "segments must be >= 1");

		let s = segments;
		let n_per_axis = s + 1;

		// ---- 混合进制 stride ----
		let mut grid_stride = [0usize; DIM];
		{
			let mut acc = 1usize;
			for d in 0..DIM {
				grid_stride[d] = acc;
				acc *= n_per_axis;
			}
		}
		let lin_of = |coords: &[usize; DIM]| -> usize {
			let mut l = 0usize;
			for d in 0..DIM {
				l += coords[d] * grid_stride[d];
			}
			l
		};

		// ---------- 1. 顶点：表面网格点投影到球面 ----------
		let total_grid = n_per_axis.pow(DIM as u32);
		let mut idx_map = std::collections::HashMap::with_capacity(total_grid);
		// let mut idx_map: Vec<usize> = Vec::with_capacity(total_grid);
		let mut vertexes: Vec<SVector<Num, DIM>> = Vec::new();

		let s_num = Num::from_usize(s).unwrap();
		let two = Num::from_usize(2).unwrap();

		for lin in 0..total_grid {
			let mut coords = [0usize; DIM];
			for d in 0..DIM {
				coords[d] = (lin / grid_stride[d]) % n_per_axis;
			}

			// 仅保留超立方体表面的点
			let on_surface = (0..DIM).any(|d| coords[d] == 0 || coords[d] == s);
			if !on_surface {
				continue;
			}

			// 网格坐标 -> [-1,1]^DIM
			let mut p = SVector::<Num, DIM>::zeros();
			for d in 0..DIM {
				let t = Num::from_usize(coords[d]).unwrap() / s_num.clone();
				p[d] = t * two.clone() - Num::one();
			}

			// 径向投影到单位球面
			let unit = p.clone() / p.norm();
			// idx_map.push(vertexes.len());
			idx_map.insert(lin, vertexes.len());
			vertexes.push(center.clone() + unit * radius.clone());
		}

		// ---------- 2. 面：排列剖分 (DIM-1)-cube ----------
		let mut faces: Vec<[usize; DIM]> = Vec::new();
		let m = DIM - 1;
		let perms = permutations(m);

		for axis in 0..DIM {
			for &fixed in &[0usize, s] {
				let rest: Vec<usize> = (0..DIM).filter(|&d| d != axis).collect();

				// 枚举 (DIM-1)-cube 里的每个 unit cell
				let n_cells = s.pow(m as u32);
				for cell in 0..n_cells {
					let mut start = [0usize; DIM];
					start[axis] = fixed;
					for k in 0..m {
						start[rest[k]] = (cell / s.pow(k as u32)) % s;
					}

					// 每个排列 => 一个 (DIM-1)-单纯形
					for perm in &perms {
						let mut coords = start;
						let mut face_arr = [0usize; DIM];

						face_arr[0] = idx_map[&lin_of(&coords)];
						for (step, &k) in perm.iter().enumerate() {
							coords[rest[k]] += 1;
							face_arr[step + 1] = idx_map[&lin_of(&coords)];
						}

						let raw_n = generalized_cross(face_arr.iter().map(|i|vertexes[*i].clone()));
						let mut fc = SVector::<Num, DIM>::zeros();
						for &vi in &face_arr {
							fc += vertexes[vi].clone();
						}
						fc /= Num::from_usize(DIM).unwrap();
						if raw_n.dot(&(fc - center.clone())) < Num::zero() {
							// face_arr.reverse();
							face_arr.swap(0, 1);
						}

						faces.push(face_arr);
					}
				}
			}
		}

		return ShapeStruct::new(vertexes, faces);
	}

	
}

fn permutations(n: usize) -> Vec<Vec<usize>> {
	let mut result = Vec::new();
	let mut cur: Vec<usize> = (0..n).collect();
	permute_rec(&mut cur, 0, &mut result);
	result
}

fn permute_rec(arr: &mut Vec<usize>, k: usize, out: &mut Vec<Vec<usize>>) {
	if k == arr.len() {
		out.push(arr.clone());
		return;
	}
	for i in k..arr.len() {
		arr.swap(k, i);
		permute_rec(arr, k + 1, out);
		arr.swap(k, i);
	}
}

#[cfg(test)]
mod tests_shape {
    use super::*;
    use nalgebra::{Vector2, Vector3, Vector4};

    const EPS: f64 = 1e-9;
    const EPS_LOOSE: f64 = 1e-6;

	fn factorial_usize(v:usize)->usize{
		factorial(v).unwrap()
	}

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < EPS
    }
    fn approx_tol(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // ============================================================
    // 1. 工具函数
    // ============================================================

    #[test]
    fn test_factorial() {
        assert_eq!(factorial_usize(0), 1);
        assert_eq!(factorial_usize(1), 1);
        assert_eq!(factorial_usize(2), 2);
        assert_eq!(factorial_usize(3), 6);
        assert_eq!(factorial_usize(4), 24);
        assert_eq!(factorial_usize(5), 120);
    }

    #[test]
    fn test_permutations() {
        assert_eq!(permutations(0), vec![Vec::<usize>::new()]);
        assert_eq!(permutations(1), vec![vec![0]]);

        let p2 = permutations(2);
        assert_eq!(p2.len(), 2);
        assert!(p2.contains(&vec![0, 1]));
        assert!(p2.contains(&vec![1, 0]));

        let p3 = permutations(3);
        assert_eq!(p3.len(), 6);
        for p in &p3 {
            let mut sorted = p.clone();
            sorted.sort();
            assert_eq!(sorted, vec![0, 1, 2]);
        }
    }
	fn normalize_or_zero<Num: RealField, const DIM: usize>(
    v: SVector<Num, DIM>,
	) -> SVector<Num, DIM> {
		let norm = v.norm();
		if norm == Num::zero() {
			v
		} else {
			v / norm
		}
	}

    #[test]
    fn test_normalize_or_zero() {
        let v = Vector3::new(3.0_f64, 0.0, 4.0);
        let u = normalize_or_zero(v);
        assert!(approx(u.norm(), 1.0));
        assert!(approx(u.x, 0.6));
        assert!(approx(u.z, 0.8));

        let z = Vector3::<f64>::zeros();
        assert_eq!(normalize_or_zero(z), z);
    }

    // ============================================================
    // 2. face_normal_raw / face_normal
    // ============================================================

	fn face_normal<Num:RealField,const DIM:usize>(vertexes: &[SVector<Num, DIM>], edge: &[usize;DIM]) -> SVector<Num, DIM> {
		normalize_or_zero(generalized_cross(edge.iter().map(|i|vertexes[*i].clone())))
	}
    #[test]
    fn test_face_normal_2d_axis_aligned() {
        let vs = [
            Vector2::new(0.0_f64, 0.0),
            Vector2::new(1.0, 0.0),
            Vector2::new(0.0, 1.0),
        ];
        // 边 (0,0)->(1,0) 法向平行 y 轴
        let n = face_normal(&vs, &[0, 1]);
        assert!(approx(n.x, 0.0));
        assert!(approx(n.y.abs(), 1.0));

        // 边 (0,0)->(0,1) 法向平行 x 轴
        let n = face_normal(&vs, &[0, 2]);
        assert!(approx(n.x.abs(), 1.0));
        assert!(approx(n.y, 0.0));
    }

    #[test]
    fn test_face_normal_2d_diagonal() {
        let vs = [Vector2::new(0.0_f64, 0.0), Vector2::new(1.0, 1.0)];
        let n = face_normal(&vs, &[0, 1]);
        let s = std::f64::consts::FRAC_1_SQRT_2;
        assert!(approx(n.x.abs(), s));
        assert!(approx(n.y.abs(), s));
        let edge = vs[1] - vs[0];
        assert!(approx(n.dot(&edge), 0.0)); // 与边垂直
    }

    #[test]
    fn test_face_normal_3d_xy_plane() {
        let vs = [
            Vector3::new(0.0_f64, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let n = face_normal(&vs, &[0, 1, 2]);
        assert!(approx(n.x, 0.0));
        assert!(approx(n.y, 0.0));
        assert!(approx(n.z.abs(), 1.0));
    }

    #[test]
    fn test_face_normal_3d_oblique_plane() {
        // 平面 z = x + y，法向沿 (1, 1, -1)
        let vs = [
            Vector3::new(0.0_f64, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
        ];
        let n = face_normal(&vs, &[0, 1, 2]);
        let expected = Vector3::new(1.0, 1.0, -1.0).normalize();
        assert!(approx(n.dot(&expected).abs(), 1.0));
    }

    #[test]
    fn test_face_normal_degenerate() {
        // 三点共线 -> 零法向
        let vs = [
            Vector3::new(0.0_f64, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
        ];
        let n = face_normal(&vs, &[0, 1, 2]);
        assert!(approx(n.norm(), 0.0));
    }

    #[test]
    fn test_face_normal_4d() {
        // 4D 中的 3-单纯形，落在超平面 x3 = 0 内
        let vs = [
            Vector4::new(0.0_f64, 0.0, 0.0, 0.0),
            Vector4::new(1.0, 0.0, 0.0, 0.0),
            Vector4::new(0.0, 1.0, 0.0, 0.0),
            Vector4::new(0.0, 0.0, 1.0, 0.0),
        ];
        let n = face_normal(&vs, &[0, 1, 2, 3]);
        assert!(approx(n.x, 0.0));
        assert!(approx(n.y, 0.0));
        assert!(approx(n.z, 0.0));
        assert!(approx(n.w.abs(), 1.0));
    }

    // ============================================================
    // 3. ShapeEdge::compute_normal / update_normal
    // ============================================================

    #[test]
    fn test_shape_edge_compute_normal_no_inside() {
        let vs = [Vector2::new(0.0_f64, 0.0), Vector2::new(1.0, 0.0)];
        // let face = ShapeEdge::<f64, 2> {
        //     vertexes_index: vec![[0, 1]],
        //     normal_vector: Vector2::zeros(),
        // };
		let face =SurfaceStruct::<f64,2>::new(&vs,[0,1]);
        let n = face.edge_vector.clone().normalize();
        assert!(approx(n.norm(), 1.0));
    }

    #[test]
    fn test_shape_face_compute_normal_with_inside() {
        // 面在 z=0，传入的 inside 在 -z 侧 => 修正后应朝 +z
        let vs = [
            Vector3::new(0.0_f64, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
		let face = SurfaceStruct::<f64,3>::new(&vs,[0,1,2]);

        let _inside = Vector3::new(0.3, 0.3, -1.0);
        let n = face.edge_vector.clone().normalize();
        assert!(n.z > 0.0, "expected +z, got {:?}", n);
    }

    #[test]
    fn test_shape_face_update_normal() {
        let vs = [Vector2::new(0.0_f64, 0.0), Vector2::new(1.0, 0.0)];
		let face = SurfaceStruct::<f64,2>::new(&vs,[0,1]);
        assert!(approx(face.edge_vector.norm(), 1.0));
    }

    // ============================================================
    // 4. ShapePoly::compute_volume / update_volume
    // ============================================================

    #[test]
    fn test_volume_2d_unit_square() {
        // let poly = ShapePoly::<f64, 2> {
        //     vertexes: vec![
        //         Vector2::new(0.0, 0.0),
        //         Vector2::new(1.0, 0.0),
        //         Vector2::new(1.0, 1.0),
        //         Vector2::new(0.0, 1.0),
        //     ],
        //     faces: vec![[0, 1], [1, 2], [2, 3], [3, 0]],
        //     volume: 0.0,
        // };
		let poly = ShapeStruct::<f64, 2>::new(
            vec![
                Vector2::new(0.0, 0.0),
                Vector2::new(1.0, 0.0),
                Vector2::new(1.0, 1.0),
                Vector2::new(0.0, 1.0),
            ],
            vec![[0, 1], [1, 2], [2, 3], [3, 0]]
        );
        assert!(approx(poly.volume, 1.0));
    }

    #[test]
    fn test_volume_2d_triangle() {
        // let poly = ShapePoly::<f64, 2> {
        //     vertexes: vec![
        //         Vector2::new(0.0, 0.0),
        //         Vector2::new(1.0, 0.0),
        //         Vector2::new(0.0, 1.0),
        //     ],
        //     faces: vec![[0, 1], [1, 2], [2, 0]],
        //     volume: 0.0,
        // };
		let poly = ShapeStruct::<f64, 2>::new(
			vec![
				Vector2::new(0.0, 0.0),
				Vector2::new(1.0, 0.0),
				Vector2::new(0.0, 1.0),
			],
			vec![[0, 1], [1, 2], [2, 0]]
		);

        assert!(approx(poly.volume, 0.5));
    }

    #[test]
    fn test_volume_3d_tetrahedron() {
        // let poly = ShapePoly::<f64, 3> {
        //     vertexes: vec![
        //         Vector3::new(0.0, 0.0, 0.0),
        //         Vector3::new(1.0, 0.0, 0.0),
        //         Vector3::new(0.0, 1.0, 0.0),
        //         Vector3::new(0.0, 0.0, 1.0),
        //     ],
        //     faces: vec![
        //         [0, 2, 1],
        //         [0, 1, 3],
        //         [0, 3, 2],
        //         [1, 2, 3],
        //     ],
        //     volume: 0.0,
        // };
		let poly = ShapeStruct::<f64, 3>::new(
			vec![
				Vector3::new(0.0, 0.0, 0.0),
				Vector3::new(1.0, 0.0, 0.0),
				Vector3::new(0.0, 1.0, 0.0),
				Vector3::new(0.0, 0.0, 1.0),
			],
			vec![
				[0, 2, 1],
				[0, 1, 3],
				[0, 3, 2],
				[1, 2, 3],
			]
		);
        assert!(approx(poly.volume, 1.0 / 6.0));
    }

    #[test]
    fn test_volume_3d_tetrahedron_translated() {
        // 平移不应改变体积
        let t = Vector3::new(5.0, -3.0, 7.0);
        // let poly = ShapePoly::<f64, 3> {
        //     vertexes: vec![
        //         Vector3::new(0.0, 0.0, 0.0) + t,
        //         Vector3::new(1.0, 0.0, 0.0) + t,
        //         Vector3::new(0.0, 1.0, 0.0) + t,
        //         Vector3::new(0.0, 0.0, 1.0) + t,
        //     ],
        //     faces: vec![
        //         [0, 2, 1],
        //         [0, 1, 3],
        //         [0, 3, 2],
        //         [1, 2, 3],
        //     ],
        //     volume: 0.0,
        // };

		let poly = ShapeStruct::<f64, 3>::new(
			vec![
				Vector3::new(0.0, 0.0, 0.0) + t,
				Vector3::new(1.0, 0.0, 0.0) + t,
				Vector3::new(0.0, 1.0, 0.0) + t,
				Vector3::new(0.0, 0.0, 1.0) + t,
			],
			vec![
				[0, 2, 1],
				[0, 1, 3],
				[0, 3, 2],
				[1, 2, 3],
			]
		);

        assert!(approx_tol(poly.volume, 1.0 / 6.0, EPS_LOOSE));
    }

    #[test]
    fn test_volume_update() {
        // let mut poly = ShapePoly::<f64, 2> {
        //     vertexes: vec![
        //         Vector2::new(0.0, 0.0),
        //         Vector2::new(2.0, 0.0),
        //         Vector2::new(2.0, 2.0),
        //         Vector2::new(0.0, 2.0),
        //     ],
        //     faces: vec![[0, 1], [1, 2], [2, 3], [3, 0]],
        //     volume: 0.0,
        // };

		let poly=ShapeStruct::<f64,2>::new(

			vec![
				Vector2::new(0.0, 0.0),
				Vector2::new(2.0, 0.0),
				Vector2::new(2.0, 2.0),
				Vector2::new(0.0, 2.0),
			],
			vec![[0, 1], [1, 2], [2, 3], [3, 0]]
		);
        assert!(approx(poly.volume, 4.0));
    }


    // ============================================================
    // 6. from_ball_grid
    // ============================================================

    #[test]
    fn test_from_ball_2d_seg1() {
        let poly = ShapeStruct::from_ball_grid(Vector2::new(0.0_f64, 0.0), 1.0, 1);
        assert_eq!(poly.vertexes.len(), 4);
        assert_eq!(poly.surfaces.len(), 4);
        for v in &poly.vertexes {
            assert!(approx(v.norm(), 1.0));
        }
        // 内接正方形，对角线长 2，面积 2
        assert!(approx(poly.volume, 2.0));
    }

    #[test]
    fn test_from_ball_2d_converges_to_pi() {
        let poly = ShapeStruct::from_ball_grid(Vector2::new(0.0_f64, 0.0), 1.0, 64);
        assert!((poly.volume - std::f64::consts::PI).abs() < 0.05);
    }

    #[test]
    fn test_from_ball_3d_seg1() {
        let poly = ShapeStruct::from_ball_grid(Vector3::new(0.0_f64, 0.0, 0.0), 1.0, 1);
        assert_eq!(poly.vertexes.len(), 8);
        assert_eq!(poly.surfaces.len(), 12);
        for v in &poly.vertexes {
            assert!(approx(v.norm(), 1.0));
        }
        // 内接立方体：边长 2/√3，体积 8/(3√3)
        let expected = 8.0 / (3.0 * 3.0_f64.sqrt());
        assert!(
            approx_tol(poly.volume, expected, 1e-9),
            "got {}, expected {}",
            poly.volume,
            expected
        );
    }

    #[test]
    fn test_from_ball_3d_converges() {
        let poly = ShapeStruct::from_ball_grid(Vector3::new(0.0_f64, 0.0, 0.0), 1.0, 8);
        let sphere_vol = 4.0 / 3.0 * std::f64::consts::PI;
        // 内接多面体体积严格小于球体积
        assert!(poly.volume < sphere_vol);
        // 8 段时通常误差 < 5%
        assert!(
            poly.volume > 0.9 * sphere_vol,
            "vol {} too small vs {}",
            poly.volume,
            sphere_vol
        );
    }

    #[test]
    fn test_from_ball_4d_seg1() {
        let poly = ShapeStruct::from_ball_grid(Vector4::new(0.0_f64, 0.0, 0.0, 0.0), 1.0, 1);
        assert_eq!(poly.vertexes.len(), 16);
        // 8 个 3-cube 面 × 3! = 6 个四面体 = 48
        assert_eq!(poly.surfaces.len(), 48);
        for v in &poly.vertexes {
            assert!(approx(v.norm(), 1.0));
        }
		for e in poly.surfaces {
			
			assert!(
				e.edge_vector.dot(&e.center)>0.0,
				"edge {:?} normal not outward",
				e.vertexes_index
			)
		}
        // 内接 4D 立方体顶点 (±0.5, ±0.5, ±0.5, ±0.5)，体积 = 1
        assert!(
            approx_tol(poly.volume, 1.0, 1e-9),
            "got {}, expected 1",
            poly.volume
        );
    }

    #[test]
    fn test_from_ball_translation_invariance() {
        let p1 = ShapeStruct::from_ball_grid(Vector3::new(0.0_f64, 0.0, 0.0), 1.0, 3);
        let p2 = ShapeStruct::from_ball_grid(Vector3::new(10.0, -5.0, 20.0), 1.0, 3);
        assert!(approx_tol(p1.volume, p2.volume, 1e-9));
    }

    // ============================================================
    // 7. 综合：法向一律朝外
    // ============================================================

    #[test]
    fn test_from_ball_all_faces_outward_2d() {
        let center = Vector2::new(0.0_f64, 0.0);
        let poly = ShapeStruct::from_ball_grid(center, 1.0, 4);
        for face in &poly.surfaces {
            let mut fc = Vector2::zeros();
            for i in face.vertexes_index {
                fc += poly.vertexes[i];
            }
            fc /= 2.0;
            let outward = fc - center;
            let raw_n = generalized_cross(face.vertexes_index.iter().map(|i|poly.vertexes[*i].clone()));
            assert!(
                raw_n.dot(&outward) > 0.0,
                "face {:?} normal not outward",
                face.vertexes_index
            );
        }
    }

    #[test]
    fn test_from_ball_all_faces_outward_3d() {
        let center = Vector3::new(0.0_f64, 0.0, 0.0);
        let poly = ShapeStruct::from_ball_grid(center, 1.0, 2);
        for face in &poly.surfaces {
            let mut fc = Vector3::zeros();
            for i in face.vertexes_index {
                fc += poly.vertexes[i];
            }
            fc /= 3.0;
            let outward = fc - center;
            let raw_n = generalized_cross(face.vertexes_index.iter().map(|i|poly.vertexes[*i].clone()));
            assert!(
                raw_n.dot(&outward) > 0.0,
                "face {:?} not outward",
                face.vertexes_index
            );
        }
    }
}
