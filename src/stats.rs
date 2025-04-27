use wacky_bag::derive_add_traits;

use crate::{num::Num, vec2_fix::Vec2Fix};


#[derive(Default,Clone,Copy,Debug)]
pub struct Mass(pub Num);

derive_add_traits!(Mass);

#[derive(Default,Clone,Copy,Debug)]
pub struct Pos(pub Vec2Fix);

derive_add_traits!(Pos);

#[derive(Default,Clone,Copy,Debug)]
pub struct Vel(pub Vec2Fix);

derive_add_traits!(Vel);


#[derive(Default,Clone,Copy,Debug)]
pub struct Dir(pub Num);

derive_add_traits!(Dir);


#[derive(Default,Clone,Copy,Debug)]
pub struct Agv(pub Num);

derive_add_traits!(Agv);


#[derive(Default,Clone,Copy,Debug)]
pub struct Energy(pub Num);

derive_add_traits!(Energy);


#[derive(Default,Clone,Copy,Debug)]
pub struct Momentum(pub Vec2Fix);

derive_add_traits!(Momentum);

#[derive(Default,Clone,Copy,Debug)]
pub struct DirVec(pub Vec2Fix);

#[derive(Default, Clone, Copy, Debug)]
pub struct Kinetic(pub Num);

derive_add_traits!(Kinetic);