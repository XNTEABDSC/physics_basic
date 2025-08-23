use wacky_bag::derive_add_traits;
use wacky_bag_fixed::vec_fix::VecFix;

use crate::{num::Num};


#[derive(Default,Clone,Copy,Debug)]
pub struct Mass(pub Num);

//derive_add_traits!{Mass}

#[derive(Default,Clone,Copy,Debug)]
pub struct Pos<const DIM:usize>(pub VecFix<DIM>);

//derive_add_traits!(Pos<1>);
//derive_add_traits!(Pos<2>);
//derive_add_traits!(Pos<3>);
//derive_add_traits!(Pos<4>);

#[derive(Default,Clone,Copy,Debug)]
pub struct Vel<const DIM:usize>(pub VecFix<DIM>);

//derive_add_traits!(Vel);


#[derive(Default,Clone,Copy,Debug)]
pub struct Dir(pub Num);

//derive_add_traits!(Dir);


#[derive(Default,Clone,Copy,Debug)]
pub struct Agv(pub Num);

//derive_add_traits!(Agv);


#[derive(Default,Clone,Copy,Debug)]
pub struct Energy(pub Num);

//derive_add_traits!(Energy);


#[derive(Default,Clone,Copy,Debug)]
pub struct Momentum(pub VecFix2);

//derive_add_traits!(Momentum);

#[derive(Default,Clone,Copy,Debug)]
pub struct DirVec(pub VecFix2);

#[derive(Default, Clone, Copy, Debug)]
pub struct Kinetic(pub Num);

//derive_add_traits!(Kinetic);

