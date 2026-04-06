
use nalgebra::{RealField, SVector, Scalar};
use num_traits::Zero;
// use wacky_bag_fixed::vec_fix::VecFix;

use derive_more::{Add,AddAssign,Sub,SubAssign,Neg};

// use crate::{num::Num};





#[derive(Default,Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign,Neg)]
pub struct Mass<Num>(pub Num);

/// Time pass of this object for this iteration of simulation
#[derive(Default,Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign,Neg)]
pub struct TimePass<Num>(pub Num);

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

impl<Num: std::ops::Neg+Scalar+simba::scalar::ClosedNeg, const DIM: usize> std::ops::Neg for $type <Num, DIM> {
    type Output=Self;

	fn neg(self) -> Self::Output {
		Self(self.0.neg())
	}
}
    };
}


//derive_add_traits!{Mass}

#[derive(Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign)]
pub struct Pos<Num,const DIM:usize>(pub SVector<Num,DIM>);
	// where Num:RealField;

derive_for_s_vector!{Pos}

//derive_add_traits!(Pos<1>);
//derive_add_traits!(Pos<2>);
//derive_add_traits!(Pos<3>);
//derive_add_traits!(Pos<4>);

#[derive(Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign)]
pub struct Vel<Num,const DIM:usize>(pub SVector<Num,DIM>);

derive_for_s_vector!{Vel}

//derive_add_traits!(Vel);


#[derive(Clone,Copy,Debug)]
pub struct DirVec<Num,const DIM:usize>(pub SVector<Num,DIM>);

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

#[derive(Default,Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign,Neg)]
pub struct Energy<Num>(pub Num);



//derive_add_traits!(Energy);


#[derive(Clone, Copy, Debug, Add, AddAssign, Sub, SubAssign)]
pub struct Momentum<Num,const DIM:usize>(pub SVector<Num,DIM>);


derive_for_s_vector!{Momentum}

//derive_add_traits!(Momentum);

#[derive(Default, Clone, Copy, Debug,Add,AddAssign,Sub,SubAssign,Neg)]
pub struct Kinetic<Num>(pub Num);

//derive_add_traits!(Kinetic);

