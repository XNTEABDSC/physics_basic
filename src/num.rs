// pub use wacky_bag_fixed::num::*;

use nalgebra::{Const, DefaultAllocator, DimMin, DimName, ToTypenum, allocator::Allocator};

use crate::rotation::DimToSoDim;


pub trait DimNameTrait:
	where Self: DimName + DimMin<Self,Output = Self> + DimToSoDim + ToTypenum,
    DefaultAllocator: Allocator< <Self as DimToSoDim>::SoDim > + Allocator<<Self as DimToSoDim>::SoDim,<Self as DimToSoDim>::SoDim>,
{

}

impl<const D:usize> DimNameTrait for Const<D> 
	where Const<D>:ToTypenum+DimToSoDim+DimMin<Self,Output = Self>,
	DefaultAllocator: Allocator< <Self as DimToSoDim>::SoDim > + Allocator<<Self as DimToSoDim>::SoDim,<Self as DimToSoDim>::SoDim>,
{
	
}