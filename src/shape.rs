use wacky_bag_fixed::num::Num;



#[derive(Default, Clone, Copy, Debug)]
pub struct ShapeCircle{
    pub radius:Num,
    pub radius_sq:Num,
}

impl ShapeCircle {
    pub fn with_radius(radius:Num)->Self {
        return Self{
            radius,
            radius_sq:radius*radius
        };
    }
}