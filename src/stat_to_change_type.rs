use std::{marker::PhantomData, ops::AddAssign};

use frunk::Func;
use num_traits::Zero;
use wacky_bag::utils::{h_list_helpers::{HTypeMapP, HZip}, type_fn::TypeFunc};

pub trait StatToChangeType<Marker>{
	type ChangeType
		// where Self:AddAssign<Self::ChangeType>,
		// 	Self::ChangeType:Default
	;
}
#[derive(Debug,Default,Clone, Copy)]
pub struct StatToChangeTypeAdd;

impl<T> StatToChangeType<StatToChangeTypeAdd> for T
where T:AddAssign<T>+Zero
{
	type ChangeType=Self;
}

pub struct MapStatToChangeTypeZ;

impl<T,M> TypeFunc<(T,M)> for MapStatToChangeTypeZ 
	where T:StatToChangeType<M>
{
	type Output=T::ChangeType;
}

impl<T,M> Func<(PhantomData<T>,PhantomData<M>)> for MapStatToChangeTypeZ 
	where T:StatToChangeType<M>
{
	type Output=PhantomData<T::ChangeType>;
	
	fn call(i: (PhantomData<T>,PhantomData<M>)) -> Self::Output {
		Default::default()
	}
	
}

impl<T,M> Func<PhantomData<(T,M)>> for MapStatToChangeTypeZ 
	where T:StatToChangeType<M>
{
	type Output=PhantomData<T::ChangeType>;
	
	fn call(i: PhantomData<(T,M)>) -> Self::Output {
		Default::default()
	}
	
}

pub type HMapStatToChangeTypeZ<T,M>=HTypeMapP<HZip<T,M>,MapStatToChangeTypeZ>;