// use wacky_bag_fixed::num::Num;


use std::array;

use frunk::{HList, Poly, hlist, hlist_pat};
use nalgebra::allocator::Allocator;
use nalgebra::{Const, DefaultAllocator, DimMin, DimName, RealField, SMatrix, SVector};
use num_traits::Zero;
use wacky_bag::utils::d_sphere_volume::d_sphere_volume_by_radius_pow;
use wacky_bag::utils::factorial::factorial;
// use wacky_bag_hlist::h_extend_by_fn::h_extend_by_fn_ref;
use wacky_bag_hlist::h_list_helpers::{HToMut, HToRef, SetMut, Sum};

use crate::stats::{Kinetic, Mass, Momentum, Pos, TimePass, Vel, Volume, mass_vel_2_kinetic};
use crate::rotation::{AllocatorDimVMSoDimVM, AngularInertia, AngularKinetic, AngularMomentum, AngularVel, ConstDimToSoDimT, DimSquare, DimSquareSo, DimToSoDim, Rotation, RotationMatrix, angular_kinetic_from_inertia_agv, angular_momentum_from_omega, angular_velocity_from_momentum};


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

#[derive(Debug)]
pub struct ShapeSurface<Num:RealField,const DIM:usize>{
	pub vertexes_index:[usize;DIM],
	pub edge_vector:SVector<Num,DIM>,
    pub normal_vector:SVector<Num,DIM>,
    pub area:Num,
	pub center:SVector<Num,DIM>
}
#[derive(Debug)]
pub struct ShapePoly<Num:RealField,const DIM:usize>{
	pub vertexes:Vec<SVector<Num,DIM>>,
	pub surfaces:Vec<ShapeSurface<Num,DIM>>,
	pub volume:Num
}

/// 计算一个面的“未归一化”法向量。
/// DIM=3 时就是三角形叉积；
/// DIM=2 时就是线段的垂直向量；
/// 一般 DIM 下使用余子式/广义叉积。
pub fn generalized_cross<Num: RealField, const DIM: usize>(
    vertexes: &[SVector<Num, DIM>],
    edge: &[usize; DIM],
) -> SVector<Num, DIM> {
    assert!(DIM >= 2);

    let p0 = vertexes[edge[0]].clone();

    // DIM 维空间中，(DIM-1)-单纯形有 DIM 个顶点，
    // 因此有 DIM-1 条边向量。
    let m = DIM - 1;
    let mut edges = Vec::with_capacity(m);

    for i in 0..m {
        edges.push(vertexes[edge[i + 1]].clone() - p0.clone());
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

impl<Num:RealField,const DIM:usize> ShapeSurface<Num,DIM> {
	pub fn new(vertexes:&[SVector<Num,DIM>],edge:[usize;DIM])->Self{
		let edge_vector=generalized_cross(vertexes,&edge);
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

impl<Num:RealField,const DIM:usize> ShapePoly<Num,DIM>
where Const<DIM>:DimSquare
{
	pub fn new(vertexes:Vec<SVector<Num,DIM>>,edges:Vec<[usize;DIM]>)->Self{
		let surfaces=edges.into_iter().map(|e|ShapeSurface::new(&vertexes,e)).collect::<Vec<_>>();
		Self{
			volume:Self::volume(&vertexes,&surfaces),
			vertexes,
			surfaces,
		}
	}
	fn volume(vertexes:&[SVector<Num,DIM>],edges:&[ShapeSurface<Num,DIM>]) -> Num 
	where Const<DIM>:DimSquare
	{
        let mut sum = Num::zero();

        for edges in edges {
            let mat = SMatrix::<Num, DIM, DIM>::from_fn(|r, c| {
                vertexes[edges.vertexes_index[r]][c].clone()
            });

            sum += mat.determinant();
        }

        let fact = Num::from_usize(factorial(DIM).unwrap()).unwrap();

        (sum / fact).abs()
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
	) -> ShapePoly<Num, DIM> {
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

						let raw_n = generalized_cross(&vertexes, &face_arr);
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

		return ShapePoly::new(vertexes, faces);
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
		normalize_or_zero(generalized_cross(vertexes, edge))
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
		let face =ShapeSurface::<f64,2>::new(&vs,[0,1]);
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
		let face = ShapeSurface::<f64,3>::new(&vs,[0,1,2]);

        let _inside = Vector3::new(0.3, 0.3, -1.0);
        let n = face.edge_vector.clone().normalize();
        assert!(n.z > 0.0, "expected +z, got {:?}", n);
    }

    #[test]
    fn test_shape_face_update_normal() {
        let vs = [Vector2::new(0.0_f64, 0.0), Vector2::new(1.0, 0.0)];
		let face = ShapeSurface::<f64,2>::new(&vs,[0,1]);
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
		let poly = ShapePoly::<f64, 2>::new(
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
		let poly = ShapePoly::<f64, 2>::new(
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
		let poly = ShapePoly::<f64, 3>::new(
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

		let poly = ShapePoly::<f64, 3>::new(
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

		let poly=ShapePoly::<f64,2>::new(

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
        let poly = ShapePoly::from_ball_grid(Vector2::new(0.0_f64, 0.0), 1.0, 1);
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
        let poly = ShapePoly::from_ball_grid(Vector2::new(0.0_f64, 0.0), 1.0, 64);
        assert!((poly.volume - std::f64::consts::PI).abs() < 0.05);
    }

    #[test]
    fn test_from_ball_3d_seg1() {
        let poly = ShapePoly::from_ball_grid(Vector3::new(0.0_f64, 0.0, 0.0), 1.0, 1);
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
        let poly = ShapePoly::from_ball_grid(Vector3::new(0.0_f64, 0.0, 0.0), 1.0, 8);
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
        let poly = ShapePoly::from_ball_grid(Vector4::new(0.0_f64, 0.0, 0.0, 0.0), 1.0, 1);
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
        let p1 = ShapePoly::from_ball_grid(Vector3::new(0.0_f64, 0.0, 0.0), 1.0, 3);
        let p2 = ShapePoly::from_ball_grid(Vector3::new(10.0, -5.0, 20.0), 1.0, 3);
        assert!(approx_tol(p1.volume, p2.volume, 1e-9));
    }

    // ============================================================
    // 7. 综合：法向一律朝外
    // ============================================================

    #[test]
    fn test_from_ball_all_faces_outward_2d() {
        let center = Vector2::new(0.0_f64, 0.0);
        let poly = ShapePoly::from_ball_grid(center, 1.0, 4);
        for face in &poly.surfaces {
            let mut fc = Vector2::zeros();
            for i in face.vertexes_index {
                fc += poly.vertexes[i];
            }
            fc /= 2.0;
            let outward = fc - center;
            let raw_n = generalized_cross(&poly.vertexes, &face.vertexes_index);
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
        let poly = ShapePoly::from_ball_grid(center, 1.0, 2);
        for face in &poly.surfaces {
            let mut fc = Vector3::zeros();
            for i in face.vertexes_index {
                fc += poly.vertexes[i];
            }
            fc /= 3.0;
            let outward = fc - center;
            let raw_n = generalized_cross(&poly.vertexes, &face.vertexes_index);
            assert!(
                raw_n.dot(&outward) > 0.0,
                "face {:?} not outward",
                face.vertexes_index
            );
        }
    }
}

/** 
`
HList!(
	TimePass<Num>,
	Mass<Num>,
	Pos<Num,DIM>,
	Momentum<Num,DIM>,
	AngularInertia<Num,DIM>,
	AngularMomentum<Num,DIM>,
	Rotation<Num,DIM>
)
`
*/


pub type PhyBodyBasic<Num,const DIM:usize>=HList!(
	TimePass<Num>,
	Mass<Num>,
	Pos<Num,DIM>,
	Momentum<Num,DIM>,
	AngularInertia<Num,DIM>,
	AngularMomentum<Num,DIM>,
	Rotation<Num,DIM>
);

/**
`
HList!(
	Mass<Num>,
	Pos<Num,DIM>,
	Momentum<Num,DIM>,
	AngularInertia<Num,DIM>,
	AngularMomentum<Num,DIM>,
	Rotation<Num,DIM>
)
`
*/

/**
```
Mass<Num>,
Pos<Num,DIM>,
Momentum<Num,DIM>,
AngularInertia<Num,DIM>,
AngularMomentum<Num,DIM>,
Rotation<Num,DIM>
```
 */
pub type PhyBodyBasicStat<Num,const DIM:usize>=HList!(
	Mass<Num>,
	Pos<Num,DIM>,
	Momentum<Num,DIM>,
	AngularInertia<Num,DIM>,
	AngularMomentum<Num,DIM>,
	Rotation<Num,DIM>
);
/**
```
TimePass<Num>,
Mass<Num>,
Pos<Num,DIM>,
Vel<Num,DIM>,
Momentum<Num,DIM>,
Kinetic<Num>,
AngularInertia<Num,DIM>,
AngularMomentum<Num,DIM>,
AngularVel<Num,DIM>,
AngularKinetic<Num>,
Rotation<Num,DIM>,	
RotationMatrix<Num,DIM>
```
 */
pub type PhyBodyFull<Num,const DIM:usize>=HList!(
	TimePass<Num>,
	Mass<Num>,
	Pos<Num,DIM>,
	Vel<Num,DIM>,
	Momentum<Num,DIM>,
	Kinetic<Num>,
	AngularInertia<Num,DIM>,
	AngularMomentum<Num,DIM>,
	AngularVel<Num,DIM>,
	AngularKinetic<Num>,
	Rotation<Num,DIM>,	
	RotationMatrix<Num,DIM>
);


pub type CalculatePositionStateInput<Num,const DIM:usize>=HList!(Mass<Num>,Momentum<Num,DIM>);

pub type CalculatePositionStateOutput<Num,const DIM:usize>=HList!(Vel<Num,DIM>,Kinetic<Num>);

/// ([`Mass`], [`Momentum`]) -> ( [`Vel`], [`Kinetic`])
pub fn calculate_position_state<'a,Num,const DIM:usize>(hlist_pat![mass,momentum] : HToRef<'a,CalculatePositionStateInput<Num,DIM>>)
->CalculatePositionStateOutput<Num,DIM>
where
	Num:RealField+Copy,
	// Const<DIM>: DimNameToSoDimName + DimName,
	// DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>+Allocator<DimNameToSoDimNameType<DIM>>,
    // DimNameToSoDimNameType<DIM>:
    //     DimMin<DimNameToSoDimNameType<DIM>, Output = DimNameToSoDimNameType<DIM>>,

{
	// let agv=angular_velocity_from_momentum(hlist![agi.clone(),&agm]).unwrap_or_else(||Zero::zero());
	let vel=Vel(momentum.0/mass.0);
	let hlist_pat![kinetic]=mass_vel_2_kinetic(hlist![mass,&vel]);
	// let ag_kinetic=angular_kinetic_from_inertia_agv(hlist![agi.clone(),&agv]);
	hlist![vel,kinetic]
}

/// [`Rotation`] -> [`RotationMatrix`]
pub fn calculate_rotation_matrix<Num:RealField+Copy,const DIM:usize>(hlist_pat![rotation]:HToRef<HList!(Rotation<Num,DIM>)>)
->HList!(RotationMatrix<Num,DIM>)
where
	Const<DIM>: DimSquare+DimToSoDim,
	DefaultAllocator:Allocator<ConstDimToSoDimT<DIM>>,
{
	hlist![RotationMatrix::from_so(rotation)]
}

/// convert [`Vel`] change to [`Momentum`]
pub fn collect_position_state_det_pos_momentum_change_vel<Num:RealField+Copy,const DIM:usize>(
	hlist_pat![mass] : HToRef<HList!(Mass<Num>)>,
	hlist_pat![vel] : HList!(Vel<Num,DIM>))
	->HList!(Momentum<Num,DIM>)
{
	hlist![Momentum(vel.0*mass.0)]
}

pub type CalculateAngularStateInput<Num,const DIM:usize>=HList!(AngularInertia<Num,DIM>,AngularMomentum<Num,DIM>);

pub type CalculateAngularStateOutput<Num,const DIM:usize>=HList!(AngularVel<Num,DIM>,AngularKinetic<Num>);

/// ([`AngularInertia`], [`AngularMomentum`]) -> ([`AngularVel`], [`AngularKinetic`])
pub fn calculate_angular_state<Num,const DIM:usize>(hlist_pat![agi,agm] : HToRef<CalculateAngularStateInput<Num,DIM>>)
->CalculateAngularStateOutput<Num,DIM>
where
	Num:RealField+Copy,
	Const<DIM>: DimSquareSo,
    ConstDimToSoDimT<DIM>: DimSquare,
    DefaultAllocator: AllocatorDimVMSoDimVM<Const<DIM>>,

{
	let hlist_pat![agv]=angular_velocity_from_momentum(hlist![agi.clone(),&agm]).unwrap_or_else(||hlist![Zero::zero()]);
	let hlist_pat![ag_kinetic]=angular_kinetic_from_inertia_agv(hlist![&agi,&agv]);
	hlist![agv,ag_kinetic]
}

/// convert [`AngularVel`] change to [`AngularMomentum`]
pub fn collect_angular_state_det_agi_agm_change_agv<Num,const DIM:usize>(
	hlist_pat![agi]:HToRef<HList!(AngularInertia<Num,DIM>)>,
	hlist_pat![agv]:HList!(AngularVel<Num,DIM>),
)->HList!(AngularMomentum<Num,DIM>)
where
	Num:RealField+Copy,
	Const<DIM>: DimToSoDim + DimName,
	DefaultAllocator: Allocator<ConstDimToSoDimT<DIM>, ConstDimToSoDimT<DIM>>+Allocator<ConstDimToSoDimT<DIM>>,
    ConstDimToSoDimT<DIM>:
        DimMin<ConstDimToSoDimT<DIM>, Output = ConstDimToSoDimT<DIM>>,

{
	angular_momentum_from_omega(hlist![agi,agv])
}


pub type CalculateBodyStateInput<Num,const DIM:usize>=Sum<CalculatePositionStateInput<Num,DIM>,CalculateAngularStateInput<Num,DIM>>;

pub type CalculateBodyStateOutput<Num,const DIM:usize>=Sum<CalculatePositionStateOutput<Num,DIM>,CalculateAngularStateOutput<Num,DIM>>;

/// [`calculate_position_state`] + [`calculate_angular_state`]
pub fn calculate_body_state<Num,const DIM:usize>(values : HToRef<CalculateBodyStateInput<Num,DIM>>)
->CalculateBodyStateOutput<Num,DIM>
where
	Num:RealField+Copy,
	Const<DIM>: DimSquareSo,
    ConstDimToSoDimT<DIM>: DimSquare,
    DefaultAllocator: AllocatorDimVMSoDimVM<Const<DIM>>,

{
	let (a,b)=values.sculpt();
	calculate_position_state(a)+calculate_angular_state(b)
	// h_extend_by_fn_ref(, calculate_rotation_matrix)
}

/// [`PhyBodyBasic`] -> [`PhyBodyFull`] by [`calculate_body_state`] + [`calculate_rotation_matrix`]
pub fn calculate_body_state_full<Num,const DIM:usize>(s:PhyBodyBasic<Num,DIM>)->PhyBodyFull<Num,DIM>
where
	Num:RealField+Copy,
	Const<DIM>: DimSquareSo,
    ConstDimToSoDimT<DIM>: DimSquare,
    DefaultAllocator: AllocatorDimVMSoDimVM<Const<DIM>>,
{
	//hlist_pat![time,mass,pos,momentum,agi,agm,rot]

	let n=calculate_body_state(s.to_ref().sculpt().0);
	let s=s+n;
	// let s=h_extend_by_fn_ref(s, calculate_body_state);
	let n=calculate_rotation_matrix(s.to_ref().sculpt().0);
	(s+n).sculpt().0
	// let a=;
	// let a=h_extend_by_fn_ref(s, calculate_body_state);
	// let a=h_extend_by_fn_ref(a, calculate_rotation_matrix);
	// a.sculpt().0
}

pub fn calculate_body_state_inplace<Num,const DIM:usize>(
	values:
	// <HToRef<CalculateBodyStateInput<Num,DIM>> as Add<HToMut<CalculateBodyStateOutput<Num,DIM>>> >::Output
	Sum<HToRef<CalculateBodyStateInput<Num,DIM>>, HToMut<CalculateBodyStateOutput<Num,DIM>>>
)
where
	Num:RealField+Copy,
	Const<DIM>: DimSquareSo,
    ConstDimToSoDimT<DIM>: DimSquare,
    DefaultAllocator: AllocatorDimVMSoDimVM<Const<DIM>>,
{
	let (inputs,outputs):(HToRef<CalculateBodyStateInput<Num,DIM>>,_)=values.sculpt();
	outputs.zip(calculate_body_state(inputs)).map(Poly(SetMut));
}

pub fn calculate_body_state_full_inplace_m<Num,const DIM:usize>(input:
	HToMut<PhyBodyFull<Num,DIM>>)
where
	Num:RealField+Copy,
	Const<DIM>: DimSquareSo,
    ConstDimToSoDimT<DIM>: DimSquare,
    DefaultAllocator: AllocatorDimVMSoDimVM<Const<DIM>>,
{
	let hlist_pat![
		_time,
		mass,
		_pos,
		vel,
		momentum,
		kinetic,
		agi,
		agm,
		agv,
		agk,
		_rot,
		_rot_mat]=input;

	calculate_body_state_inplace(hlist![mass,momentum,agi,agm,vel,kinetic,agv,agk]);
}


#[cfg(test)]
mod tests{
	
    use frunk::hlist::HZippable;
    use wacky_bag_hlist::{h_list_helpers::{HMapP, MapToPhantom}};

    use crate::{body::PhyBodyBasicStat, stat_to_change_type::MapStatToChangeTypeZ};

	type Num=f32;
	const DIM:usize=3;
	#[test]
	fn test(){
		// type Changes= HMapP<<HMapP<PhyBodyBasicStat<Num,DIM>,MapToPhantom> as HZippable<_>>::Zipped,Poly<MapStatToChangeTypeZ>> ;
		// let dwa:HMapP<
		// 	<HMapP<PhyBodyBasicStat<Num,DIM>,MapToPhantom> 
		// 	as HZippable<
		// 		HMapP<HList!(StatToChangeTypeAdd,StatToChangeTypeAdd,StatToChangeTypeAdd,StatToChangeTypeAdd,StatToChangeTypeAdd,StatToChangeTypeAdd,RotationToRotationDelta),MapToPhantom>
		// 	>>::Zipped,MapStatToChangeTypeZ
		// >=Default::default();

		let _dwa2:HMapP<
			<HMapP<
				PhyBodyBasicStat<Num,DIM>
				// HList!(Mass<Num>,Pos<Num,DIM>,Momentum<Num,DIM>,Energy<Num>,)
				,MapToPhantom> 
			as HZippable<
				// HMapP<HList!(StatToChangeTypeAdd),MapToPhantom>
				_
			>>::Zipped,MapStatToChangeTypeZ
		>=Default::default();
		// type Changes=;
	}
}