// use wacky_bag_fixed::num::Num;

use nalgebra::RealField;



#[derive(Default, Clone, Copy, Debug)]
pub struct ShapeCircle<Num:RealField>{
    pub radius:Num,
    pub radius_sq:Num,
}

impl<Num:RealField+Copy+num_traits::Num> ShapeCircle<Num> {
    pub fn from_radius(radius:Num)->Self {
        return Self{
            radius,
            radius_sq:radius*radius
        };
    }
}