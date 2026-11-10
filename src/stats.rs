
use frunk::{HList, hlist, hlist_pat};
use nalgebra::{Const, DefaultAllocator, OVector, RealField, SVector, Scalar, allocator::Allocator};
use num_traits::Zero;
// use wacky_bag_fixed::vec_fix::VecFix;
use wacky_bag::utils::num_extend::NumExtends;
use derive_more::{Add,AddAssign,Sub,SubAssign,Neg};

use crate::derive_zero;

#[derive(Clone, Copy, Debug, Add, AddAssign, Sub, SubAssign, Neg, PartialEq, Eq, PartialOrd, Ord)]
pub struct Mass<Num>(pub Num);

derive_zero!(Mass);

/// Time pass of this object for this iteration of simulation
#[derive(Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign,Neg, PartialEq, Eq, PartialOrd, Ord)]
pub struct TimePass<Num>(pub Num);

derive_zero!(TimePass);

macro_rules! derive_for_s_vector {
    ($type:tt) => {
impl<Num,const DIM:usize> Default for $type <Num,DIM>
	where Num:RealField
{
    fn default() -> Self {
        Self(SVector::zeros())
    }
}

impl<Num,const DIM:usize> Zero for $type <Num,DIM> 
	where Num:RealField
{
	fn zero() -> Self {
		Self(SVector::zero())
	}

	fn is_zero(&self) -> bool {
		self.0.is_zero()
	}
}

// impl<Num: RealField, const DIM: usize> std::ops::Neg for $type <Num, DIM> {
//     type Output=Self;

// 	fn neg(self) -> Self::Output {
// 		Self(self.0.neg())
// 	}
// }
    };
}


//derive_add_traits!{Mass}

#[derive(Clone, Debug, Add, AddAssign, Sub, SubAssign, Neg, PartialEq, Eq)]
pub struct Pos<Num,const DIM:usize>(pub OVector<Num,Const<DIM>>)
where Num:RealField,DefaultAllocator:Allocator<Const<DIM>>;

impl<Num: RealField + Copy, const DIM: usize> Copy for Pos<Num, DIM>
where DefaultAllocator:Allocator<Const<DIM>,Buffer<Num>:Copy>
{
}
	// where Num:RealField;

derive_for_s_vector!{Pos}

//derive_add_traits!(Pos<1>);
//derive_add_traits!(Pos<2>);
//derive_add_traits!(Pos<3>);
//derive_add_traits!(Pos<4>);

#[derive(Clone, Debug, Add, AddAssign, Sub, SubAssign, Neg, PartialEq, Eq)]
pub struct Vel<Num,const DIM:usize>(pub OVector<Num,Const<DIM>>)
where Num:RealField,DefaultAllocator:Allocator<Const<DIM>>;

impl<Num: Copy, const DIM: usize> Copy for Vel<Num, DIM>
where Num:RealField,DefaultAllocator:Allocator<Const<DIM>,Buffer<Num>:Copy>
{
}

derive_for_s_vector!{Vel}

//derive_add_traits!(Vel);


#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DirVec<Num,const DIM:usize>(pub OVector<Num,Const<DIM>>)
where Num:RealField,DefaultAllocator:Allocator<Const<DIM>>;

impl<Num: Copy, const DIM: usize> Copy for DirVec<Num, DIM>
where Num:RealField,DefaultAllocator:Allocator<Const<DIM>,Buffer<Num>:Copy>
{
}

impl<Num,const DIM:usize> Default for DirVec<Num,DIM>
	where Num:RealField
{
    fn default() -> Self {
        Self(SVector::zeros())
    }
}

//derive_add_traits!(Dir);

/*
#[derive(Default,Clone,Copy,Debug)]
pub struct Agv<const DIM:usize>(pub CMatrix<Num,DIM,DIM>);
//derive_add_traits!(Agv);

 */

#[derive(Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign,Neg, PartialEq, Eq, PartialOrd)]
pub struct Energy<Num>(pub Num);

derive_zero!(Energy);


//derive_add_traits!(Energy);


#[derive(Clone, Debug, Add, AddAssign, Sub, SubAssign,Neg, PartialEq, Eq)]
pub struct Momentum<Num,const DIM:usize>(pub OVector<Num,Const<DIM>>)
where Num:RealField,DefaultAllocator:Allocator<Const<DIM>>;

impl<Num: Copy, const DIM: usize> Copy for Momentum<Num, DIM>
where Num:RealField,DefaultAllocator:Allocator<Const<DIM>,Buffer<Num>:Copy>
{
}


derive_for_s_vector!{Momentum}

//derive_add_traits!(Momentum);

#[derive(Clone, Copy, Debug,Add,AddAssign,Sub,SubAssign,Neg, PartialEq, Eq, PartialOrd)]
pub struct Kinetic<Num>(pub Num);

derive_zero!(Kinetic);

pub fn mass_vel_2_momentum<Num:RealField+Copy,const DIM:usize>(hlist_pat![mass,vel]:HList!(&Mass<Num>,&Vel<Num,DIM>))->HList!(Momentum<Num,DIM>) {
	
	hlist![
		// Kinetic(mass.0*velsq/Num::frac_1_2())
		Momentum(vel.0*mass.0)
	]
}

pub fn mass_vel_2_kinetic<Num:RealField+Copy,const DIM:usize>(hlist_pat![mass,vel]:HList!(&Mass<Num>,&Vel<Num,DIM>))->HList!(Kinetic<Num>) {
	let velsq=vel.0.dot(&vel.0);
	hlist![
		Kinetic(mass.0*velsq/Num::frac_1_2())
	]
}

#[derive(Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign,Neg, PartialEq, Eq, PartialOrd)]
pub struct Volume<Num>(pub Num);

derive_zero!(Volume);

#[derive(Default, Clone, Copy, Debug, PartialEq, Eq, PartialOrd)]
pub struct Density<Num>(pub Num);

//derive_add_traits!(Kinetic);

