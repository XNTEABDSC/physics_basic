use cordic::*;
use crate::num::Num;
//pub struct Vec2F(pub Num,pub Num);

pub type Vec2Fix=wacky_bag::structures::vec2::Vec2<Num>;//super::vec2::Vec2<Num>;

#[inline]
pub fn dist_sq(a: Vec2Fix)->Num {
    a*a
}

#[inline]
pub fn dist(a: Vec2Fix)->Num {
    cordic::sqrt(a*a)
}

#[inline]
pub fn dir(a: Vec2Fix)->Num {
    cordic::atan2(a[1], a[0])
}

#[inline]
pub fn from_rot(rot:Num)->Vec2Fix{
    Vec2Fix::new(cos(rot),sin(rot))
}
#[inline]
pub fn from_polar(rot:Num,dist:Num)->Vec2Fix {
    Vec2Fix::new(cos(rot)*dist,sin(rot)*dist)
}