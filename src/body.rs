// use wacky_bag_fixed::num::Num;

use frunk::{HList, hlist, hlist_pat};
use nalgebra::allocator::Allocator;
use nalgebra::{Const, DefaultAllocator, DimMin, DimName, Matrix, RealField, SMatrix};
use wacky_bag::utils::h_list_helpers::{HToMut, HToRef};

use crate::stats::{Mass,Momentum,Energy,Pos,TimePass,Vel};
use crate::rotation::{AngularInertia, AngularMomentum, AngularVel, DimNameToSoDimName, DimNameToSoDimNameType, Rotation, RotationDelta, angular_velocity_from_momentum};


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

/// `
/// HList!(
/// 	TimePass<Num>,
/// 	Mass<Num>,
/// 	Pos<Num,DIM>,
/// 	Momentum<Num,DIM>,
/// 	AngularInertia<Num,DIM>,
/// 	AngularMomentum<Num,DIM>,
/// 	Rotation<Num,DIM>
/// )
/// `
pub type PhyBodyBasic<Num,const DIM:usize>=HList!(
	TimePass<Num>,
	Mass<Num>,
	Pos<Num,DIM>,
	Momentum<Num,DIM>,
	// Energy<Num>,
	AngularInertia<Num,DIM>,
	AngularMomentum<Num,DIM>,
	Rotation<Num,DIM>
);

/// `
/// HList!(
/// 	Mass<Num>,
/// 	Pos<Num,DIM>,
/// 	Momentum<Num,DIM>,
/// 	AngularInertia<Num,DIM>,
/// 	AngularMomentum<Num,DIM>,
/// 	Rotation<Num,DIM>
/// )
/// `
pub type PhyBodyBasicStat<Num,const DIM:usize>=HList!(
	Mass<Num>,
	Pos<Num,DIM>,
	Momentum<Num,DIM>,
	// Energy<Num>,
	AngularInertia<Num,DIM>,
	AngularMomentum<Num,DIM>,
	Rotation<Num,DIM>
);

pub type PhyBodyFull<Num,const DIM:usize>=HList!(
	TimePass<Num>,
	Mass<Num>,
	Pos<Num,DIM>,
	Vel<Num,DIM>,
	Momentum<Num,DIM>,
	// Energy<Num>,
	AngularInertia<Num,DIM>,
	AngularMomentum<Num,DIM>,
	AngularVel<Num,DIM>,
	Rotation<Num,DIM>,	
);

pub type CalculateBodyStateInput<Num,const DIM:usize>=HList!(Mass<Num>,Momentum<Num,DIM>,AngularInertia<Num,DIM>,AngularMomentum<Num,DIM>);
pub type CalculateBodyStateOutput<Num,const DIM:usize>=HList!(Vel<Num,DIM>,AngularVel<Num,DIM>);

pub fn calculate_body_state<Num,const DIM:usize>(hlist_pat![mass,momentum,agi,agm]:HToRef<CalculateBodyStateInput<Num,DIM>>)
->CalculateBodyStateOutput<Num,DIM>
where
	Num:RealField+Copy,
	Const<DIM>: DimNameToSoDimName + DimName,
	DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>+Allocator<DimNameToSoDimNameType<DIM>>,
    DimNameToSoDimNameType<DIM>:
        DimMin<DimNameToSoDimNameType<DIM>, Output = DimNameToSoDimNameType<DIM>>,

{
	let agv=angular_velocity_from_momentum(hlist![agi.clone(),&agm]).unwrap_or_else(||SMatrix::<Num,DIM,DIM>::zeros());
	hlist![Vel(momentum.0/mass.0),AngularVel(agv)]
}

pub fn calculate_body_state_full<Num,const DIM:usize>(s:PhyBodyBasic<Num,DIM>)->PhyBodyFull<Num,DIM>
where
	Num:RealField+Copy,
	Const<DIM>: DimNameToSoDimName + DimName,
	DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>+Allocator<DimNameToSoDimNameType<DIM>>,
    DimNameToSoDimNameType<DIM>:
        DimMin<DimNameToSoDimNameType<DIM>, Output = DimNameToSoDimNameType<DIM>>,
{
	//hlist_pat![time,mass,pos,momentum,agi,agm,rot]

	let n=calculate_body_state(s.to_ref().sculpt().0);
	(s+n).sculpt().0
	// let agv=angular_velocity_from_momentum(hlist![agi.clone(),&agm]).unwrap_or_else(||SMatrix::<Num,DIM,DIM>::zeros());
	// hlist![
	// 	time,
	// 	mass,
	// 	pos,
	// 	Vel(momentum.0/mass.0),
	// 	momentum,
	// 	agi,
	// 	agm,
	// 	AngularVel(agv),
	// 	rot
	// ]
}

pub fn calculate_body_state_inplace<Num,const DIM:usize>(hlist_pat![
	mass,momentum,agi,agm,
	vel,agv
	]:
	HList!(&Mass<Num>,&Momentum<Num,DIM>,&AngularInertia<Num,DIM>,&AngularMomentum<Num,DIM>
	,&mut Vel<Num,DIM>,&mut AngularVel<Num,DIM>))
where
	Num:RealField+Copy,
	Const<DIM>: DimNameToSoDimName + DimName,
	DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>+Allocator<DimNameToSoDimNameType<DIM>>,
    DimNameToSoDimNameType<DIM>:
        DimMin<DimNameToSoDimNameType<DIM>, Output = DimNameToSoDimNameType<DIM>>,
{
	vel.0=momentum.0/mass.0;
	agv.0=angular_velocity_from_momentum(hlist![agi.clone(),&agm]).unwrap_or_else(||SMatrix::<Num,DIM,DIM>::zeros());
}






pub fn calculate_body_state_full_inplace_m<Num,const DIM:usize>(hlist_pat![
		time,
		mass,
		pos,
		vel,
		momentum,
		agi,
		agm,
		agv,
		rot]:
	HToMut<PhyBodyFull<Num,DIM>>)
where
	Num:RealField+Copy,
	Const<DIM>: DimNameToSoDimName + DimName,
	DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>+Allocator<DimNameToSoDimNameType<DIM>>,
    DimNameToSoDimNameType<DIM>:
        DimMin<DimNameToSoDimNameType<DIM>, Output = DimNameToSoDimNameType<DIM>>,
{
	calculate_body_state_inplace(hlist![mass,momentum,agi,agm,vel,agv]);
}










#[cfg(test)]
mod tests{
	use super::*;
    use frunk::{HList, Poly, hlist::HZippable};
    use wacky_bag::utils::{h_list_helpers::{HMapP, MapToPhantom}};

    use crate::{body::PhyBodyBasicStat, rotation::RotationToRotationDelta, stat_to_change_type::{MapStatToChangeTypeZ, StatToChangeTypeAdd}, stats::{Mass, Pos}};

	type Num=f32;
	const DIM:usize=3;
	#[test]
	fn test(){
		// type Changes= HMapP<<HMapP<PhyBodyBasicStat<Num,DIM>,MapToPhantom> as HZippable<_>>::Zipped,Poly<MapStatToChangeTypeZ>> ;
		// let dwa:HMapP<
		// 	<HMapP<PhyBodyBasicStat<Num,DIM>,MapToPhantom> 
		// 	as HZippable<
		// 		HMapP<HList!(StatToChangeTypeAdd,StatToChangeTypeAdd,StatToChangeTypeAdd,StatToChangeTypeAdd,StatToChangeTypeAdd,StatToChangeTypeAdd,RotationToRotationDelta),MapToPhantom>
		// 	>>::Zipped,MapStatToChangeTypeZ
		// >=Default::default();

		let dwa2:HMapP<
			<HMapP<
				PhyBodyBasicStat<Num,DIM>
				// HList!(Mass<Num>,Pos<Num,DIM>,Momentum<Num,DIM>,Energy<Num>,)
				,MapToPhantom> 
			as HZippable<
				// HMapP<HList!(StatToChangeTypeAdd),MapToPhantom>
				_
			>>::Zipped,MapStatToChangeTypeZ
		>=Default::default();
		// type Changes=;
	}
}