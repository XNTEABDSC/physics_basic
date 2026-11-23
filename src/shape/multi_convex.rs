use nalgebra::RealField;

use crate::shape::convex_polyhedron::ConvexPolyhedron;

pub struct MultipleConvex<Num: RealField + Copy, const DIM: usize>{
    convexes:Vec<ConvexPolyhedron<Num,DIM>>
}

