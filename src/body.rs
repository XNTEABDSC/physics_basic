// use wacky_bag_fixed::num::Num;

use frunk::{HList, hlist, hlist_pat};
use nalgebra::allocator::Allocator;
use nalgebra::{Const, DefaultAllocator, DimMin, DimName, Matrix, RealField, SMatrix};
use num_traits::Zero;
use wacky_bag::utils::h_list_helpers::{HToMut, HToRef};

use crate::stats::{Energy, Kinetic, Mass, Momentum, Pos, TimePass, Vel, kinetic_from_mass_vel};
use crate::rotation::{AngularInertia, AngularKinetic, AngularMomentum, AngularVel, DimNameToSoDimName, DimNameToSoDimNameType, Rotation, RotationDelta, angular_kinetic_from_inertia_agv, angular_velocity_from_momentum};


#[derive(Default, Clone, Copy, Debug)]
pub struct ShapeCircle<Num:RealField,const DIM:usize>{
    pub radius:Num,
    pub radius_pow:Num,
}

impl<Num:RealField+Copy+num_traits::Num,const DIM:usize> ShapeCircle<Num,DIM> {
    pub fn from_radius(radius:Num)->Self {
        return Self{
            radius,
            radius_pow:radius.powi(DIM as i32)
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
	Kinetic<Num>,
	// Energy<Num>,
	AngularInertia<Num,DIM>,
	AngularMomentum<Num,DIM>,
	AngularVel<Num,DIM>,
	AngularKinetic<Num>,
	Rotation<Num,DIM>,	
);

pub type CalculateBodyStateInput<Num,const DIM:usize>=HList!(Mass<Num>,Momentum<Num,DIM>,AngularInertia<Num,DIM>,AngularMomentum<Num,DIM>);
pub type CalculateBodyStateOutput<Num,const DIM:usize>=HList!(Vel<Num,DIM>,AngularVel<Num,DIM>,Kinetic<Num>,AngularKinetic<Num>);
pub fn calculate_body_state<Num,const DIM:usize>(hlist_pat![mass,momentum,agi,agm] : HToRef<CalculateBodyStateInput<Num,DIM>>)
->CalculateBodyStateOutput<Num,DIM>
where
	Num:RealField+Copy,
	Const<DIM>: DimNameToSoDimName + DimName,
	DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>+Allocator<DimNameToSoDimNameType<DIM>>,
    DimNameToSoDimNameType<DIM>:
        DimMin<DimNameToSoDimNameType<DIM>, Output = DimNameToSoDimNameType<DIM>>,

{
	let agv=angular_velocity_from_momentum(hlist![agi.clone(),&agm]).unwrap_or_else(||Zero::zero());
	let vel=Vel(momentum.0/mass.0);
	let kinetic=kinetic_from_mass_vel(hlist![mass,&vel]);
	let ag_kinetic=angular_kinetic_from_inertia_agv(hlist![agi.clone(),&agv]);
	hlist![vel,agv,kinetic,ag_kinetic]
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
	*agv=angular_velocity_from_momentum(hlist![agi.clone(),&agm]).unwrap_or_else(||Zero::zero());
}






pub fn calculate_body_state_full_inplace_m<Num,const DIM:usize>(hlist_pat![
		time,
		mass,
		pos,
		vel,
		momentum,
		kinetic,
		agi,
		agm,
		agv,
		agk,
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