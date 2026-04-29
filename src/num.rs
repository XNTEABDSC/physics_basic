// pub use wacky_bag_fixed::num::*;

use nalgebra::{Const, DefaultAllocator, DimMin, DimName, ToTypenum, allocator::Allocator};

use crate::rotation::DimNameToSoDimName;


pub trait DimNameTrait:
	where Self: DimName + DimMin<Self,Output = Self> + DimNameToSoDimName + ToTypenum,
    DefaultAllocator: Allocator< <Self as DimNameToSoDimName>::SoDimName > + Allocator<<Self as DimNameToSoDimName>::SoDimName,<Self as DimNameToSoDimName>::SoDimName>,
{

}

impl<const D:usize> DimNameTrait for Const<D> 
	where Const<D>:ToTypenum+DimNameToSoDimName+DimMin<Self,Output = Self>,
	DefaultAllocator: Allocator< <Self as DimNameToSoDimName>::SoDimName > + Allocator<<Self as DimNameToSoDimName>::SoDimName,<Self as DimNameToSoDimName>::SoDimName>,
{
	
}