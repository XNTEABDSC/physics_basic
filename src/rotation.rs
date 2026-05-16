use std::ops::{Add, AddAssign, Div, Mul, Rem, Sub};

use crate::{
    num::DimNameTrait, stat_to_change_type::StatToChangeType, stats::{Mass, Pos}
};
use derive_more::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use frunk::{HList, hlist, hlist_pat};
use nalgebra::{
    Const, DefaultAllocator, DimMin, DimName, LU, OMatrix, OVector, RealField, SMatrix, SVector,
    ToConst, ToTypenum, allocator::Allocator,
};
use name_type_for_fn::name_type;
use num_traits::Zero;
use typenum::{U0, U1, U2};
use wacky_bag::utils::num_extend::NumExtends;
// use simba::scalar::FixedI32F32;

pub trait DimNameToSoDimName
where
    Self: ToTypenum<
        Typenum: Sub<U1>
                     + Mul<
            <Self::Typenum as Sub<U1>>::Output,
            Output: Rem<U2, Output = U0> + Div<U2, Output: ToConst>,
        >,
    >,
{
    type SoDimName: DimName;
    const SO_DIM: usize;
}

// type ToSoSizeTExp<N>= <Prod<N,Sum<N,N1>> as PartialDiv<P2>>::Output;

impl<const DIM: usize, TNum, SubTNum1, TNum2> DimNameToSoDimName for Const<DIM>
where
    Self: ToTypenum<Typenum = TNum>,
    TNum: Sub<U1, Output = SubTNum1>
        + Mul<SubTNum1, Output: Div<U2, Output = TNum2> + Rem<U2, Output = U0>>,
    TNum2: ToConst,
{
    type SoDimName = TNum2::Const;

    const SO_DIM: usize = TNum2::Const::DIM;
}

pub type DimNameToSoDimNameType<const DIM: usize> = <Const<DIM> as DimNameToSoDimName>::SoDimName;

#[derive(Clone, Copy, Debug, Add, AddAssign, Sub, SubAssign, Neg)]
pub struct AngularVel<Num: RealField, const DIM: usize>(pub SMatrix<Num, DIM, DIM>);

impl<Num: RealField, const DIM: usize> Default for AngularVel<Num, DIM> {
	fn default() -> Self {
		Self(SMatrix::<Num, DIM, DIM>::zeros())
	}
}

impl<Num: RealField, const DIM: usize> Zero for AngularVel<Num, DIM> {
	
	fn zero() -> Self {
		Self(SMatrix::<Num, DIM, DIM>::zeros())
	}
	
	fn is_zero(&self) -> bool {
		self.0.is_zero()
	}
}

#[derive(Clone, Copy, Debug, Add, AddAssign, Sub, SubAssign, Neg)]
pub struct AngularMomentum<Num: RealField, const DIM: usize>(pub SMatrix<Num, DIM, DIM>);

impl<Num: RealField, const DIM: usize> Default for AngularMomentum<Num, DIM> {
	fn default() -> Self {
		Self(SMatrix::<Num, DIM, DIM>::zeros())
	}
}

impl<Num: RealField, const DIM: usize> Zero for AngularMomentum<Num, DIM> {
	
	fn zero() -> Self {
		Self(SMatrix::<Num, DIM, DIM>::zeros())
	}
	
	fn is_zero(&self) -> bool {
		self.0.is_zero()
	}
}

#[derive(Clone, Copy, Debug, Mul, MulAssign)]
pub struct Rotation<Num: RealField, const DIM: usize>(pub SMatrix<Num, DIM, DIM>);

impl<Num: RealField, const DIM: usize> Default for Rotation<Num, DIM> {
    fn default() -> Self {
        Self(SMatrix::<Num, DIM, DIM>::identity())
    }
}

// pub trait AllocatorForSo<const DIM:usize> = Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>,Buffer :Send+Sync>;

#[derive(Clone, Debug, Add, AddAssign, Sub, SubAssign, Neg)]
pub struct AngularInertia<Num: RealField+Copy, const DIM: usize>(
    pub OMatrix<Num, DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>,
)
where
    Const<DIM>: DimNameToSoDimName,
    DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>;

impl<Num: RealField + Copy, const DIM: usize> Default for AngularInertia<Num, DIM>
where
	Const<DIM>: DimNameToSoDimName,
	DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>
{
	fn default() -> Self {
		Self(OMatrix::<Num, DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>::zeros())
	}
}

impl<Num: RealField + Copy, const DIM: usize> Zero for AngularInertia<Num, DIM>
where
	Const<DIM>: DimNameToSoDimName,
	DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>
{
	fn zero() -> Self {
		Self(OMatrix::<Num, DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>::zeros())
	}

	fn is_zero(&self) -> bool {
		self.0.is_zero()
	}
	// fn default() -> Self {
	// 	Self(OMatrix::<Num, DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>::zeros())
	// }
}

impl<Num: RealField+Copy, const DIM: usize> Copy for AngularInertia<Num,DIM>
where
    Const<DIM>: DimNameToSoDimName,
    DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>,Buffer<Num> : Copy>
{
	
}


#[derive(Clone, Copy, Debug, Add, AddAssign, Sub, SubAssign, Neg)]
pub struct RotationDelta<Num: RealField, const DIM: usize>(pub SMatrix<Num, DIM, DIM>);

impl<Num: RealField, const DIM: usize> AddAssign<RotationDelta<Num, DIM>> for Rotation<Num, DIM>
where
    Const<DIM>: DimMin<Const<DIM>, Output = Const<DIM>>,
{
    fn add_assign(&mut self, rhs: RotationDelta<Num, DIM>) {
        *self *= rhs.0.exp();
    }
}

pub struct RotationToRotationDelta;

impl<Num: RealField, const DIM: usize> StatToChangeType<RotationToRotationDelta> for Rotation<Num,DIM> {
	type ChangeType=RotationDelta<Num,DIM>;
}

impl<Num: RealField, const DIM: usize> AddAssign<AngularVel<Num, DIM>> for Rotation<Num, DIM>
where
    Const<DIM>: DimMin<Const<DIM>, Output = Const<DIM>>,
{
    fn add_assign(&mut self, rhs: AngularVel<Num, DIM>) {
        *self *= rhs.0.exp();
    }
}

impl<Num: RealField, const DIM: usize> Default for RotationDelta<Num, DIM> {
    fn default() -> Self {
        Self(SMatrix::<Num, DIM, DIM>::zeros())
    }
}

impl<Num: RealField, const DIM: usize> Zero for RotationDelta<Num, DIM> {
	
	fn zero() -> Self {
		Self(SMatrix::<Num, DIM, DIM>::zeros())
	}
	
	fn is_zero(&self) -> bool {
		self.0.is_zero()
	}
}


pub fn angular_vel_to_rotation<Num: RealField, const DIM: usize>(
    agv: &AngularVel<Num, DIM>,
    dt: Num,
) -> SMatrix<Num, DIM, DIM>
where
    Const<DIM>: DimMin<Const<DIM>, Output = Const<DIM>>,
{
    (&agv.0 * dt).exp()
}

pub fn so_to_vec<T: RealField + Copy, const DIM: usize>(
    omega: &SMatrix<T, DIM, DIM>,
) -> OVector<T, DimNameToSoDimNameType<DIM>>
where
    Const<DIM>: DimNameToSoDimName,
    DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>>,
{
    // const N: usize = D * (D - 1) / 2;
    let mut vec = OVector::<T, DimNameToSoDimNameType<DIM>>::zeros();
    let mut idx = 0;
    for i in 0..DIM {
        for j in i + 1..DIM {
            vec[idx] = omega[(i, j)];
            idx += 1;
        }
    }
    vec
}
/// 将so(n)基下的向量转换回反对称矩阵
pub fn vec_to_so<T: RealField + Copy, const DIM: usize>(
    vec: &OVector<T, DimNameToSoDimNameType<DIM>>,
) -> SMatrix<T, DIM, DIM>
where
    Const<DIM>: DimNameToSoDimName,
    DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>>
        + Allocator<Const<DIM>, Const<DIM>, Buffer<T> = nalgebra::ArrayStorage<T, DIM, DIM>>,
{
    let mut mat = SMatrix::<T, DIM, DIM>::zeros();
    let mut idx = 0;
    for i in 0..DIM {
        for j in i + 1..DIM {
            let val = vec[idx];
            mat[(i, j)] = val;
            mat[(j, i)] = -val;
            idx += 1;
        }
    }
    mat
}

/// 计算单个质点的轨道角动量贡献矩阵（在so(n)基下）
/// 公式: M_{ab,ij} = m * (δ_{b,i} r_a r_j - δ_{b,j} r_a r_i - δ_{a,i} r_b r_j + δ_{a,j} r_b r_i)
fn orbit_matrix<T: RealField + Copy, const DIM: usize>(
    mass: T,
    r: &SVector<T, DIM>,
) -> OMatrix<T, DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>
where
    Const<DIM>: DimNameToSoDimName + DimName,
    DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>>
        + Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>,
{
    let mut m = OMatrix::<T, DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>::zeros();
    let mut out_idx = 0;
    for a in 0..DIM {
        for b in a + 1..DIM {
            let mut in_idx = 0;
            for i in 0..DIM {
                for j in i + 1..DIM {
                    let term1 = if b == i {
                        mass * r[a] * r[j]
                    } else {
                        T::zero()
                    };
                    let term2 = if b == j {
                        -mass * r[a] * r[i]
                    } else {
                        T::zero()
                    };
                    let term3 = if a == i {
                        -mass * r[b] * r[j]
                    } else {
                        T::zero()
                    };
                    let term4 = if a == j {
                        mass * r[b] * r[i]
                    } else {
                        T::zero()
                    };
                    m[(out_idx, in_idx)] = term1 + term2 + term3 + term4;
                    in_idx += 1;
                }
            }
            out_idx += 1;
        }
    }
    m
}

/// 计算总转动惯量矩阵（在so(n)基下）
pub fn total_inertia_matrix<'a, T: RealField + Copy, const DIM: usize>(
    // objects: &[(T, SVector<T, D>, SMatrix<T, {to_so_size!(D)}, {to_so_size!(D)}>)],
    objects: impl IntoIterator<Item = HList!(&'a Mass<T>, &'a Pos<T, DIM>, &'a AngularInertia<T, DIM>)>,
) -> AngularInertia<T, DIM>
where
    Const<DIM>: DimNameToSoDimName + DimName,
    DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>>
        + Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>,
{
    let res = objects.into_iter().fold(
        OMatrix::<T, DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>::zeros(),
        move |acc, hlist_pat![mass, pos, inertia]| {
            // acc + orbit_matrix(mass.0, &pos.0) + inertia.0
            let a = orbit_matrix(mass.0, &pos.0);
            let b = &inertia.0;
            let i = a + b;
            <OMatrix<T, DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>> as Add<_>>::add(
                acc, &i,
            )
        },
    );
    return AngularInertia(res);
}

/// 根据角速度计算总角动量（反对称矩阵形式）
pub fn angular_momentum_from_omega<T: RealField + Copy, const DIM: usize>(
    hlist_pat![total_inertia, omega]: HList!(&AngularInertia<T,DIM>,&AngularVel<T, DIM>),
) -> HList!(AngularMomentum<T, DIM>)
//SMatrix<T, DIM, DIM>
where
    Const<DIM>: DimNameToSoDimName + DimName,
    DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>>
        + Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>,
{
    let omega_vec = so_to_vec(&omega.0);
    let total_inertia = &total_inertia.0;
    let l_vec = total_inertia * omega_vec;
    hlist![AngularMomentum(vec_to_so(&l_vec))]
}

/// 根据总角动量计算角速度（反对称矩阵形式）
/// 返回 `None` 若转动惯量矩阵奇异（刚体可绕某些轴自由旋转）
pub fn angular_velocity_from_momentum_<T: RealField + Copy, D2>(
    angular_momentum_so: &OVector<T, D2>,
    total_inertia: OMatrix<T, D2, D2>,
) -> Option<OVector<T, D2>>
where
    D2: DimMin<D2, Output = D2>,
    DefaultAllocator: Allocator<D2> + Allocator<D2, D2>,
{
    // let angular_momentum=angular_momentum.0;
    // let angular_momentum_so = so_to_vec(&angular_momentum);
    let decomp = LU::new(total_inertia);
    if decomp.is_invertible() {
        let omega_vec = decomp.solve(angular_momentum_so)?;
        Some(omega_vec)
    } else {
        None
    }
}

// pub where awd=Const<DIM>: DimNameToSoDimName + DimName;
// pub trait DimSO = DimNameToSoDimName + DimName ;

/// 根据总角动量计算角速度（反对称矩阵形式）
/// 返回 `None` 若转动惯量矩阵奇异（刚体可绕某些轴自由旋转）
pub fn angular_velocity_from_momentum<T: RealField + Copy, const DIM: usize>(
    hlist_pat![total_inertia, angular_momentum]: HList!(AngularInertia<T,DIM>,&AngularMomentum<T, DIM>),
) -> Option<HList!(AngularVel<T,DIM>)>
where
    Const<DIM>: DimNameToSoDimName + DimName,
    DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>>
        + Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>,
    DimNameToSoDimNameType<DIM>:
        DimMin<DimNameToSoDimNameType<DIM>, Output = DimNameToSoDimNameType<DIM>>,
{
	
    angular_velocity_from_momentum_(&so_to_vec(&angular_momentum.0), total_inertia.0)
        .map(|v| hlist![AngularVel(vec_to_so(&v))])
}

/// 计算 D 维均匀超球体（半径为 r）的转动惯量。
/// 结果在 so(n) 基下是一个对角矩阵，所有对角元相等，值为 2 M r^2 / (n+2)。
pub fn sphere_inertia<T: RealField + Copy, const DIM: usize>(
    mass: T,
    radius: T,
) -> AngularInertia<T, DIM>
where
    Const<DIM>: DimNameToSoDimName + DimName,
    DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>> + Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>,
{
    let n = T::from_usize(DIM).unwrap();
    let two = T::one() + T::one();
    let factor = two / (n + two);          // 2/(n+2)
    let i0 = mass * radius * radius * factor;

    let n_so = DimNameToSoDimNameType::<DIM>::dim();
    let mut mat = OMatrix::<T, DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>::zeros();
    for i in 0..n_so {
        mat[(i, i)] = i0;
    }
    AngularInertia(mat)
}

/// 计算 D 维均匀超长方体（边平行于坐标轴，中心在原点）的转动惯量。
/// side_lengths 是一个长度为 D 的向量，给出每个维度的边长。
/// 结果在 so(n) 基下是对角矩阵，对应平面 (i,j) 的对角元为 M (l_i^2 + l_j^2) / 12。
pub fn cuboid_inertia<T: RealField + Copy, const DIM: usize>(
    mass: T,
    side_lengths: &SVector<T, DIM>,
) -> AngularInertia<T, DIM>
where
    Const<DIM>: DimNameToSoDimName + DimName,
    DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>> + Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>,
{
    let twelve = T::from_usize(12).unwrap();
    let n_so = DimNameToSoDimNameType::<DIM>::dim();
    let mut mat = OMatrix::<T, DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>::zeros();
    let mut idx = 0;
    for i in 0..DIM {
        for j in i + 1..DIM {
            let l_i = side_lengths[i];
            let l_j = side_lengths[j];
            let i_val = mass * (l_i * l_i + l_j * l_j) / twelve;
            mat[(idx, idx)] = i_val;
            idx += 1;
        }
    }
    AngularInertia(mat)
}

#[derive(Default, Clone, Copy, Debug,Add,AddAssign,Sub,SubAssign,Neg)]
pub struct AngularKinetic<Num>(pub Num);

impl <Num: Zero>Zero for AngularKinetic<Num>{
    fn zero() -> Self {
        Self(Num::zero())
    }
    fn is_zero(&self) -> bool {
        Num::is_zero(&self.0)
    }
}

/// 基于角速度向量（反对称矩阵）和转动惯量计算转动动能。
/// 动能 = 1/2 * ω_vecᵀ · I · ω_vec
pub fn angular_kinetic_from_inertia_agv<Num: RealField + Copy, const DIM: usize>(
	hlist_pat![inertia, angular_vel]: HList!(&AngularInertia<Num,DIM>,&AngularVel<Num, DIM>),
) -> HList!(AngularKinetic<Num>)
where
    Const<DIM>: DimNameToSoDimName + DimName,
    DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>> + Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>,
{
    let omega_vec = so_to_vec(&angular_vel.0);
    let i_omega = &inertia.0 * omega_vec.clone();
    let half = Num::p2();
    hlist![AngularKinetic( half * omega_vec.dot(&i_omega) )]
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use frunk::hlist;
    use nalgebra::{Const, OMatrix, SMatrix, SVector};

    // 辅助函数：创建零惯量（质点）
    fn zero_inertia<T: RealField + Copy, const D: usize>() -> AngularInertia<T, D>
    where
        Const<D>: DimNameToSoDimName + DimName,
        DefaultAllocator: Allocator<DimNameToSoDimNameType<D>>
            + Allocator<DimNameToSoDimNameType<D>, DimNameToSoDimNameType<D>>,
    {
        AngularInertia(OMatrix::<
            T,
            DimNameToSoDimNameType<D>,
            DimNameToSoDimNameType<D>,
        >::zeros())
    }

    #[test]
    fn test_2d_single_particle() {
        type T = f64;
        const DIM: usize = 2;

        let mass = Mass(1.0);
        let pos = Pos(SVector::<T, DIM>::new(1.0, 0.0));
        let inertia = zero_inertia::<T, DIM>();

        // 角速度：绕垂直轴旋转，ω = 1
        let mut omega_mat = SMatrix::<T, DIM, DIM>::zeros();
        omega_mat[(0, 1)] = -1.0; // ω_{12} = -ω (so_to_vec 取上三角)
        omega_mat[(1, 0)] = 1.0;
        let angular_vel = AngularVel(omega_mat);

        // 计算总惯量矩阵
        let objects = vec![hlist!(&mass, &pos, &inertia)];
        let total_inertia = total_inertia_matrix(objects);

        // 计算角动量
        let hlist_pat![angular_momentum] = angular_momentum_from_omega(hlist!(&total_inertia, &angular_vel));

        // 理论值：对于 2D 单质点，角动量应为 m * r^2 * ω = 1 * 1^2 * 1 = 1
        // 反对称矩阵只有一个独立分量，即 L_{12} = -1？需要根据约定：
        // 角速度 ω 对应矩阵 [0, -ω; ω, 0]，角动量应为 [0, -L; L, 0]，其中 L = Iω = 1
        // 所以 L_{12} = -1
        let expected_l12 = 1.0;
        assert_relative_eq!(angular_momentum.0[(0, 1)], expected_l12, epsilon = 1e-12);
        assert_relative_eq!(angular_momentum.0[(1, 0)], -expected_l12, epsilon = 1e-12);

        // 反解角速度
        let hlist_pat![AngularVel(solved_omega)] = angular_velocity_from_momentum(hlist!(total_inertia, &angular_momentum))
            .expect("Should be invertible");
        assert_relative_eq!(solved_omega[(0, 1)], angular_vel.0[(0, 1)], epsilon = 1e-12);
        assert_relative_eq!(solved_omega[(1, 0)], angular_vel.0[(1, 0)], epsilon = 1e-12);
    }

    #[test]
    fn test_3d_single_particle() {
        type T = f64;
        const DIM: usize = 3;

        let mass = Mass(1.0);
        let pos = Pos(SVector::<T, DIM>::new(1.0, 0.0, 0.0));
        let inertia = zero_inertia::<T, DIM>();

        let mut omega_mat = SMatrix::<T, DIM, DIM>::zeros();
        omega_mat[(0, 1)] = -1.0;
        omega_mat[(1, 0)] = 1.0;
        let angular_vel = AngularVel(omega_mat);

        let objects = vec![hlist!(&mass, &pos, &inertia)];
        let total_inertia = total_inertia_matrix(objects);
        let hlist_pat![angular_momentum] = angular_momentum_from_omega(hlist!(&total_inertia, &angular_vel));

        // 单质点系统奇异，无法唯一确定角速度，应返回 None
        let result = angular_velocity_from_momentum(hlist!(total_inertia, &angular_momentum));
        assert!(result.is_none());
    }

    #[test]
    fn test_3d_two_particles() {
        type T = f64;
        const DIM: usize = 3;

        let mass1 = Mass(1.0);
        let pos1 = Pos(SVector::<T, DIM>::new(1.0, 0.0, 0.0));
        let mass2 = Mass(2.0);
        let pos2 = Pos(SVector::<T, DIM>::new(0.0, 1.0, 0.0));
        let inertia1 = zero_inertia::<T, DIM>();
        let inertia2 = zero_inertia::<T, DIM>();

        let objects = vec![
            hlist!(&mass1, &pos1, &inertia1),
            hlist!(&mass2, &pos2, &inertia2),
        ];
        let total_inertia = total_inertia_matrix(objects);

        // 角速度：绕 z 轴旋转 ω_z = 1
        let mut omega_mat = SMatrix::<T, DIM, DIM>::zeros();
        omega_mat[(0, 1)] = -1.0;
        omega_mat[(1, 0)] = 1.0;
        let angular_vel = AngularVel(omega_mat);

        let hlist_pat![angular_momentum] = angular_momentum_from_omega(hlist!(&total_inertia, &angular_vel));

        // 理论角动量：L = m1 r1^2 ω + m2 r2^2 ω = (1*1^2 + 2*1^2) * 1 = 3，方向沿 z
        // 因此反对称矩阵应为 3 * omega_mat
        let expected = -3.0 * omega_mat;
        assert_relative_eq!(angular_momentum.0, expected, epsilon = 1e-12);

        // 反解
        let hlist_pat![AngularVel(solved_omega)] = angular_velocity_from_momentum(hlist!(total_inertia, &angular_momentum))
            .expect("Should be invertible");
        assert_relative_eq!(solved_omega, omega_mat, epsilon = 1e-12);
    }

    #[test]
    fn test_3d_with_self_inertia() {
        type T = f64;
        const DIM: usize = 3;

        let mass = Mass(1.0);
        let pos = Pos(SVector::<T, DIM>::new(0.0, 0.0, 0.0)); // 质心在原点
        // 自转惯量（在 so(3) 基下）：假设绕 z 轴转动惯量为 2，其它轴为 1
        // so(3) 基向量顺序： (0,1), (0,2), (1,2) 对应角速度分量 ω_z, ω_y, ω_x? 需要根据代码约定。
        // 实际轨道惯量矩阵与自转惯量直接相加，测试仅验证线性关系，我们使用对角矩阵
        let mut inertia_mat =
            OMatrix::<T, DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>::zeros();
        // 假设对角元素为 [I_xy, I_xz, I_yz] 分别对应 ω_z, ω_y, ω_x 的系数
        inertia_mat[(0, 0)] = 2.0; // I_zz?
        inertia_mat[(1, 1)] = 1.0;
        inertia_mat[(2, 2)] = 1.0;
        let inertia = AngularInertia(inertia_mat);

        let objects = vec![hlist!(&mass, &pos, &inertia)];
        let total_inertia = total_inertia_matrix(objects);

        // 角速度：绕 z 轴 ω_z = 1
        let mut omega_mat = SMatrix::<T, DIM, DIM>::zeros();
        omega_mat[(0, 1)] = -1.0;
        omega_mat[(1, 0)] = 1.0;
        let angular_vel = AngularVel(omega_mat);

        let hlist_pat![angular_momentum] = angular_momentum_from_omega(hlist!(&total_inertia, &angular_vel));

        // 预期角动量：自转贡献 L_self = I * ω_vec，其中 ω_vec = so_to_vec(omega) = [-1, 0, 0]（因为 (0,1) 元素为 -1）
        // 所以 L_vec = [2*(-1), 1*0, 1*0] = [-2, 0, 0]，对应反对称矩阵：只有 (0,1) 分量为 -2
        let expected_omega_mat = 2.0 * omega_mat; // 因为 ω_mat 的 (0,1) 是 -1，乘以 2 得 -2
        assert_relative_eq!(angular_momentum.0, expected_omega_mat, epsilon = 1e-12);

        // 反解
        let hlist_pat![AngularVel(solved_omega)] = angular_velocity_from_momentum(hlist!(total_inertia, &angular_momentum))
            .expect("Should be invertible");
        assert_relative_eq!(solved_omega, omega_mat, epsilon = 1e-12);
    }

    #[test]
    fn test_singular_inertia() {
        type T = f64;
        const DIM: usize = 3;

        // 所有质点共线，导致惯性张量奇异
        let mass1 = Mass(1.0);
        let pos1 = Pos(SVector::<T, DIM>::new(1.0, 0.0, 0.0));
        let mass2 = Mass(1.0);
        let pos2 = Pos(SVector::<T, DIM>::new(2.0, 0.0, 0.0));
        let inertia1 = zero_inertia::<T, DIM>();
        let inertia2 = zero_inertia::<T, DIM>();

        let objects = vec![
            hlist!(&mass1, &pos1, &inertia1),
            hlist!(&mass2, &pos2, &inertia2),
        ];
        let total_inertia = total_inertia_matrix(objects);

        // 任意角动量（例如绕 z 轴）
        let mut l_mat = SMatrix::<T, DIM, DIM>::zeros();
        l_mat[(0, 1)] = -1.0;
        l_mat[(1, 0)] = 1.0;
        let angular_momentum = AngularMomentum(l_mat);

        // 应该无法解出角速度（惯性矩阵奇异）
        let result = angular_velocity_from_momentum(hlist!(total_inertia, &angular_momentum));
        assert!(result.is_none());
    }

	

    type T = f64;
    const D3: usize = 3;
    const D4: usize = 4;

    #[test]
    fn test_sphere_3d() {
        let mass = 1.0;
        let radius = 1.0;
        let inertia = sphere_inertia::<T, D3>(mass, radius);
        // 3D 球体: I = 2/5 M R^2 = 0.4
        let expected = 0.4;
        // so(3) 基有 3 个对角元，应全部相等
        for i in 0..3 {
            assert_relative_eq!(inertia.0[(i, i)], expected, epsilon = 1e-12);
        }
    }

    #[test]
    fn test_cuboid_3d() {
        let mass = 1.0;
        let sides = SVector::<T, D3>::new(2.0, 2.0, 2.0); // 立方体边长 2
        let inertia = cuboid_inertia::<T, D3>(mass, &sides);
        // 对于立方体，绕 x、y、z 轴的转动惯量均为 M*(a^2+a^2)/12 = 2M a^2/12 = M a^2/6
        // a=2 => 4/6=2/3≈0.6667，但注意 so(3) 基顺序：(0,1) 对应绕 z 轴，(0,2) 绕 y 轴，(1,2) 绕 x 轴
        let expected = 2.0 / 3.0; // 0.666666...
        for i in 0..3 {
            assert_relative_eq!(inertia.0[(i, i)], expected, epsilon = 1e-12);
        }
    }

    #[test]
    fn test_cuboid_3d_rect() {
        let mass = 1.0;
        let sides = SVector::<T, D3>::new(2.0, 4.0, 6.0);
        let inertia = cuboid_inertia::<T, D3>(mass, &sides);
        // 手动计算各平面惯量
        let i01 = mass * (2.0*2.0 + 4.0*4.0) / 12.0; // (2^2+4^2)/12 = (4+16)/12=20/12=1.6667
        let i02 = mass * (2.0*2.0 + 6.0*6.0) / 12.0; // (4+36)/12=40/12=3.3333
        let i12 = mass * (4.0*4.0 + 6.0*6.0) / 12.0; // (16+36)/12=52/12=4.3333
        assert_relative_eq!(inertia.0[(0, 0)], i01, epsilon = 1e-12);
        assert_relative_eq!(inertia.0[(1, 1)], i02, epsilon = 1e-12);
        assert_relative_eq!(inertia.0[(2, 2)], i12, epsilon = 1e-12);
    }
}
