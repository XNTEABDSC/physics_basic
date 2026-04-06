use std::ops::Add;
use std::ops::Mul;

use frunk::HList;
use frunk::hlist_pat;
use nalgebra::Const;
use nalgebra::{Const, DimMin, LU, RealField, SMatrix, SVector, ToConst, ToTypenum};
use derive_more::{Add,AddAssign,Sub,SubAssign,Neg};
use typenum::{Min, N1, P2, P3, P4, PartialDiv, Prod, Sum};
use crate::stats::Pos;

use crate::stats::Mass;
// use simba::scalar::FixedI32F32;

macro_rules! to_so_size {
	($d:ident) => {
		$d*($d-1)/2
	};
}

#[derive(Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign,Neg)]
pub struct AngularVel<Num:RealField,const DIM:usize>(pub SMatrix<Num,DIM,DIM>);

#[derive(Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign,Neg)]
pub struct AngularMomentum<Num:RealField,const DIM:usize>(pub SMatrix<Num,DIM,DIM>);

#[derive(Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign,Neg)]
pub struct Rotation<Num:RealField,const DIM:usize>(pub SMatrix<Num,DIM,DIM>);

pub struct AngularInertia<Num:RealField,const DIM_TEST:usize>(pub SMatrix<Num,{DIM_TEST*(DIM_TEST-1)/2},{DIM_TEST*(DIM_TEST-1)/2}>) where [(); DIM_TEST*(DIM_TEST-1)/2]:;

pub fn angular_vel_to_rotation<Num:RealField,const DIM: usize>(
    agv: &AngularVel<Num,DIM>,
    dt: Num,
) -> SMatrix<Num,DIM,DIM> 
	where Const<DIM>: DimMin<Const<DIM>,Output = Const<DIM>>
{
	(&agv.0 * dt) .exp()
}

// use nalgebra::{SMatrix, SVector, OMatrix, Matrix, Const, Dim, DefaultAllocator, allocator::Allocator};
// use nalgebra::linalg::LU;
// use num_traits::Num;
// use simba::scalar::RealField;
// use std::ops::Mul;

// fn so_to_vec<T: RealField, const D: usize>(omega: &SMatrix<T, D, D>) -> SVector<T, { D*(D-1)/2 }> {
//     let n = D*(D-1)/2;
//     let mut vec = SVector::<T, { D*(D-1)/2 }>::zeros();
//     let mut idx = 0;
//     for i in 0..D {
//         for j in i+1..D {
//             vec[idx] = omega[(i, j)]; // 上三角元素
//             idx += 1;
//         }
//     }
//     vec
// }

// type ToSoSize<N>=<Prod<N,Sum<N,N1>> as PartialDiv<P2>>::Output;

// type Test1=ToSoSize<P4>;
// #[test]
// fn test(){
// 	// let a:Test1=Default::default();
// 	println!("{}",<Test1 as typenum::Integer>::to_i32());
// 	// let v1:Sum<P3,N1>=Default::default();
// }

// fn so_to_vec<T: RealField, D>(omega: &SMatrix<T, { <D as ToConst>::Const::DIM }, D>) -> SVector<T, <ToSoSize<D> as ToConst>::Const>
// where
// 	D:ToConst<>
// {
//     let n = D::dim();          // 维度 D 的 usize 值
//     let mut vec = SVector::<T, SoSize<D>>::zeros();
//     let mut idx = 0;
//     for i in 0..n {
//         for j in i+1..n {
//             vec[idx] = omega[(i, j)];
//             idx += 1;
//         }
//     }
//     vec
// }


// pub const N<const D:usize>:usize=D*D;


pub fn so_to_vec<T: RealField+Copy, const D: usize>(omega: &SMatrix<T, D, D>) -> SVector<T, { D*(D-1)/2 }> {
    // const N: usize = D * (D - 1) / 2;
    let mut vec = SVector::<T, { D*(D-1)/2 }>::zeros();
    let mut idx = 0;
    for i in 0..D {
        for j in i+1..D {
            vec[idx] = omega[(i, j)];
            idx += 1;
        }
    }
    vec
}

/// 将so(n)基下的向量转换回反对称矩阵
pub fn vec_to_so<T: RealField+Copy, const D: usize>(vec: &SVector<T, { D*(D-1)/2 }>) -> SMatrix<T, D, D> {
    // const N: usize = D * (D - 1) / 2;
    let mut mat = SMatrix::<T, D, D>::zeros();
    let mut idx = 0;
    for i in 0..D {
        for j in i+1..D {
            let val = vec[idx];
            mat[(i, j)] = val;
            mat[(j, i)] = -val;
            idx += 1;
        }
    }
    mat
}

/// 计算单个质点的轨道角动量贡献矩阵（在so(n)基下）
fn orbit_matrix<T: RealField+Copy, const D: usize>(mass: T, r: &SVector<T, D>) -> SMatrix<T, {to_so_size!(D)}, {to_so_size!(D)}> {
    // const N: usize = D * (D - 1) / 2;
    let mut m = SMatrix::<T, {to_so_size!(D)}, {to_so_size!(D)}>::zeros();
    let mut out_idx = 0;
    for a in 0..D {
        for b in a+1..D {
            let mut in_idx = 0;
            for i in 0..D {
                for j in i+1..D {
                    // 轨道角动量对基的偏导
                    let term1 = if b == i {
                        mass * r[a] * r[j]
                    } else if b == j {
                        -mass * r[a] * r[i]
                    } else {
                        T::zero()
                    };
                    let term2 = if a == i {
                        -mass * r[b] * r[j]
                    } else if a == j {
                        mass * r[b] * r[i]
                    } else {
                        T::zero()
                    };
                    m[(out_idx, in_idx)] = term1 + term2;
                    in_idx += 1;
                }
            }
            out_idx += 1;
        }
    }
    m
}

/// 计算总转动惯量矩阵（在so(n)基下）
fn total_inertia_matrix<'a,T: RealField+Copy, const D: usize>(
    // objects: &[(T, SVector<T, D>, SMatrix<T, {to_so_size!(D)}, {to_so_size!(D)}>)],
	objects:impl IntoIterator<Item = HList!(&'a Mass<T>,&'a Pos<T,D>,&'a AngularInertia<T,D>)>
) -> SMatrix<T, {to_so_size!(D)}, {to_so_size!(D)}> {
    // const N: usize = {D * (D - 1) / 2};
    // let mut total = SMatrix::<T, {to_so_size!(D)}, {to_so_size!(D)}>::zeros();
    // for hlist_pat!(mass, pos, inertia) in objects {
    //     total += orbit_matrix(mass.0, &pos.0) + inertia.0;
    // }
    // total
	objects.into_iter().fold(SMatrix::<T, {to_so_size!(D)}, {to_so_size!(D)}>::zeros(), move |acc,hlist_pat![mass,pos,inertia]|{
		// acc + orbit_matrix(mass.0, &pos.0) + inertia.0
		let i=orbit_matrix(mass.0, &pos.0) + inertia.0;
		<SMatrix<T, {to_so_size!(D)}, {to_so_size!(D)}> as Add<_>>::add(acc, &i)
	})
}

/// 根据角速度计算总角动量（反对称矩阵形式）
pub fn angular_momentum_from_omega<T: RealField+Copy, const D: usize>(
    // objects: impl IntoIterator<Item = HList!(&'a Mass<T>,&'a Pos<T,D>,&'a AngularInteria<T,D>)>,
	// total_inertia:&AngularInertia<T,D>,
    // omega: &SMatrix<T, D, D>,
	hlist_pat![total_inertia,omega]:HList!(&AngularInertia<T,D>,&SMatrix<T, D, D>)
) -> SMatrix<T, D, D>
	where [(); D*(D-1)/2]:
{
    let omega_vec = so_to_vec(omega);
	let total_inertia=&total_inertia.0;
    // let total_inertia = total_inertia_matrix(objects);
    let l_vec = total_inertia*omega_vec;
		// <SMatrix<T, {to_so_size!(D)}, {to_so_size!(D)}> as Mul<&SVector<T, {to_so_size!(D)}>>>::mul(total_inertia, &omega_vec);
	//total_inertia * omega_vec;
    vec_to_so(&l_vec)
}


/// 根据总角动量计算角速度（反对称矩阵形式）
/// 返回 `None` 若转动惯量矩阵奇异（刚体可绕某些轴自由旋转）
pub fn angular_velocity_from_momentum_<T: RealField+Copy, const D2: usize>(
    // objects: &[(T, SVector<T, D>, SMatrix<T, { D*(D-1)/2 }, { D*(D-1)/2 }>)],
    // angular_momentum: &SMatrix<T, D, D>,

	
	// hlist_pat![total_inertia,angular_momentum]:HList!(&AngularInertia<T,D2>,&AngularMomentum<T, D>)
	angular_momentum_so:SVector<T,D2>,
	total_inertia:SMatrix<T,D2,D2>
) -> Option<SVector<T, D2>>
	where Const<D2>:DimMin<Const<D2>,Output = Const<D2>>,
		// nalgebra::DefaultAllocator: 
		// 	nalgebra::allocator::Allocator<nalgebra::Const<D2>, nalgebra::Const<D2>> 
		// 	+ nalgebra::allocator::Allocator<nalgebra::DimMinimum<nalgebra::Const<D2>, nalgebra::Const<D2>>>
	// nalgebra::Const<{D*(D-1)/2}>:DimMin<nalgebra::Const<{ D*(D-1)/2 }>>,
	// where nalgebra::Const<{ D*(D-1)/2 }>: DimMin<nalgebra::Const<{ D*(D-1)/2 }>>,
	// where nalgebra::Const<{ D*(D-1)/2 }>: ToTypenum<Typenum :Min >
{
	// let angular_momentum=angular_momentum.0;
    // let angular_momentum_so = so_to_vec(&angular_momentum);
    let decomp = LU::new(total_inertia);
    if decomp.is_invertible() {
        let omega_vec = decomp.solve(&angular_momentum_so)?;
        Some(omega_vec)
    } else {
        None
    }
}

/// 根据总角动量计算角速度（反对称矩阵形式）
/// 返回 `None` 若转动惯量矩阵奇异（刚体可绕某些轴自由旋转）
pub fn angular_velocity_from_momentum<T: RealField+Copy, const D: usize>(
    // objects: &[(T, SVector<T, D>, SMatrix<T, { D*(D-1)/2 }, { D*(D-1)/2 }>)],
    // angular_momentum: &SMatrix<T, D, D>,

	
	hlist_pat![total_inertia,angular_momentum]:HList!(&AngularInertia<T,D>,&AngularMomentum<T, D>)
) -> Option<SMatrix<T, D, D>>
	where [(); D*(D-1)/2]:,
		// nalgebra::Const<{ D*(D-1)/2 }>: DimMin<nalgebra::Const<{ D*(D-1)/2 }>>
	// nalgebra::Const<{D*(D-1)/2}>:DimMin<nalgebra::Const<{ D*(D-1)/2 }>>,
	// where nalgebra::Const<{ D*(D-1)/2 }>: DimMin<nalgebra::Const<{ D*(D-1)/2 }>>,
	// where nalgebra::Const<{ D*(D-1)/2 }>: ToTypenum<Typenum :Min >
{

	angular_velocity_from_momentum_(
		so_to_vec(&angular_momentum.0),
		total_inertia.0
	).map(|v|vec_to_so(&v))
	// let angular_momentum=angular_momentum.0;
    // let l_vec = so_to_vec(&angular_momentum);
    // let total_inertia = total_inertia.0;
    // let decomp = LU::new(total_inertia);
    // if decomp.is_invertible() {
    //     let omega_vec = decomp.solve(&l_vec)?;
    //     Some(vec_to_so(&omega_vec))
    // } else {
    //     None
    // }
}

