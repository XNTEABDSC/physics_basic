use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use wacky_bag_fixed::vec_fix::VecFix;

use derive_more::{Add,AddAssign,Sub,SubAssign,Neg};

use crate::{num::Num};





#[derive(Default,Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign,Neg)]
pub struct Mass(pub Num);
#[derive(Default,Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign,Neg)]
pub struct TimePass(pub Num);

macro_rules! derive_default_as_zeros {
    ($type:tt) => {
impl<const DIM:usize> Default for $type <DIM>{
    fn default() -> Self {
        Self(VecFix::<DIM>::zeros())
    }
}
    };
}


//derive_add_traits!{Mass}

#[derive(Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign,Neg)]
pub struct Pos<const DIM:usize>(pub VecFix<DIM>);

derive_default_as_zeros!{Pos}


//derive_add_traits!(Pos<1>);
//derive_add_traits!(Pos<2>);
//derive_add_traits!(Pos<3>);
//derive_add_traits!(Pos<4>);

#[derive(Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign,Neg)]
pub struct Vel<const DIM:usize>(pub VecFix<DIM>);

derive_default_as_zeros!{Vel}

//derive_add_traits!(Vel);


#[derive(Clone,Copy,Debug)]
pub struct DirVec<const DIM:usize>(pub VecFix<DIM>);

//derive_add_traits!(Dir);

/*
#[derive(Default,Clone,Copy,Debug)]
pub struct Agv<const DIM:usize>(pub CMatrix<Num,DIM,DIM>);
//derive_add_traits!(Agv);

 */

#[derive(Default,Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign,Neg)]
pub struct Energy(pub Num);



//derive_add_traits!(Energy);


#[derive(Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign,Neg)]
pub struct Momentum<const DIM:usize>(pub VecFix<DIM>);


derive_default_as_zeros!{Momentum}

//derive_add_traits!(Momentum);

#[derive(Default, Clone, Copy, Debug,Add,AddAssign,Sub,SubAssign,Neg)]
pub struct Kinetic(pub Num);

//derive_add_traits!(Kinetic);

