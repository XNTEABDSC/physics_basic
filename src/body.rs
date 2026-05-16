// use wacky_bag_fixed::num::Num;

use std::ops::Add;

use frunk::{HList, Poly, hlist, hlist_pat};
use nalgebra::allocator::Allocator;
use nalgebra::{Const, DefaultAllocator, DimMin, DimName, Matrix, RealField, SMatrix};
use num_traits::Zero;
use wacky_bag::utils::d_sphere_volume::d_sphere_volume_by_radius_pow;
use wacky_bag::utils::h_list_helpers::{HToMut, HToRef, SetMut, Sum};

use crate::stats::{Energy, Kinetic, Mass, Momentum, Pos, TimePass, Vel, Volume, mass_vel_2_kinetic};
use crate::rotation::{AngularInertia, AngularKinetic, AngularMomentum, AngularVel, DimNameToSoDimName, DimNameToSoDimNameType, Rotation, RotationDelta, angular_kinetic_from_inertia_agv, angular_velocity_from_momentum};


#[derive(Default, Clone, Copy, Debug)]
pub struct ShapeSphere<Num:RealField,const DIM:usize>{
    pub radius:Num,
    pub radius_pow:Num,
	pub volume:Num
}

impl<Num:RealField+Copy+num_traits::Num,const DIM:usize> ShapeSphere<Num,DIM> {
    pub fn from_radius(radius:Num)->Self {
		let radius_pow=radius.powi(DIM as i32);
        return Self{
            radius,
            radius_pow,
			volume:d_sphere_volume_by_radius_pow(radius_pow, DIM)
        };
    }

	pub fn volume(&self)->Volume<Num>{Volume(self.volume)}
}

/** 
`
HList!(
	TimePass<Num>,
	Mass<Num>,
	Pos<Num,DIM>,
	Momentum<Num,DIM>,
	AngularInertia<Num,DIM>,
	AngularMomentum<Num,DIM>,
	Rotation<Num,DIM>
)
`
*/


pub type PhyBodyBasic<Num,const DIM:usize>=HList!(
	TimePass<Num>,
	Mass<Num>,
	Pos<Num,DIM>,
	Momentum<Num,DIM>,
	AngularInertia<Num,DIM>,
	AngularMomentum<Num,DIM>,
	Rotation<Num,DIM>
);

/**
`
HList!(
	Mass<Num>,
	Pos<Num,DIM>,
	Momentum<Num,DIM>,
	AngularInertia<Num,DIM>,
	AngularMomentum<Num,DIM>,
	Rotation<Num,DIM>
)
`
*/
pub type PhyBodyBasicStat<Num,const DIM:usize>=HList!(
	Mass<Num>,
	Pos<Num,DIM>,
	Momentum<Num,DIM>,
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
	AngularInertia<Num,DIM>,
	AngularMomentum<Num,DIM>,
	AngularVel<Num,DIM>,
	AngularKinetic<Num>,
	Rotation<Num,DIM>,	
);


pub type CalculatePositionStateInput<Num,const DIM:usize>=HList!(Mass<Num>,Momentum<Num,DIM>);

pub type CalculatePositionStateOutput<Num,const DIM:usize>=HList!(Vel<Num,DIM>,Kinetic<Num>);

pub fn calculate_position_state<'a,Num,const DIM:usize>(hlist_pat![mass,momentum] : HToRef<'a,CalculatePositionStateInput<Num,DIM>>)
->CalculatePositionStateOutput<Num,DIM>
where
	Num:RealField+Copy,
	// Const<DIM>: DimNameToSoDimName + DimName,
	// DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>+Allocator<DimNameToSoDimNameType<DIM>>,
    // DimNameToSoDimNameType<DIM>:
    //     DimMin<DimNameToSoDimNameType<DIM>, Output = DimNameToSoDimNameType<DIM>>,

{
	// let agv=angular_velocity_from_momentum(hlist![agi.clone(),&agm]).unwrap_or_else(||Zero::zero());
	let vel=Vel(momentum.0/mass.0);
	let hlist_pat![kinetic]=mass_vel_2_kinetic(hlist![mass,&vel]);
	// let ag_kinetic=angular_kinetic_from_inertia_agv(hlist![agi.clone(),&agv]);
	hlist![vel,kinetic]
}


pub type CalculateAngularStateInput<Num,const DIM:usize>=HList!(AngularInertia<Num,DIM>,AngularMomentum<Num,DIM>);

pub type CalculateAngularStateOutput<Num,const DIM:usize>=HList!(AngularVel<Num,DIM>,AngularKinetic<Num>);

pub fn calculate_angular_state<Num,const DIM:usize>(hlist_pat![agi,agm] : HToRef<CalculateAngularStateInput<Num,DIM>>)
->CalculateAngularStateOutput<Num,DIM>
where
	Num:RealField+Copy,
	Const<DIM>: DimNameToSoDimName + DimName,
	DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>+Allocator<DimNameToSoDimNameType<DIM>>,
    DimNameToSoDimNameType<DIM>:
        DimMin<DimNameToSoDimNameType<DIM>, Output = DimNameToSoDimNameType<DIM>>,

{
	let hlist_pat![agv]=angular_velocity_from_momentum(hlist![agi.clone(),&agm]).unwrap_or_else(||hlist![Zero::zero()]);
	let hlist_pat![ag_kinetic]=angular_kinetic_from_inertia_agv(hlist![&agi,&agv]);
	hlist![agv,ag_kinetic]
}

pub type CalculateBodyStateInput<Num,const DIM:usize>=Sum<CalculatePositionStateInput<Num,DIM>,CalculateAngularStateInput<Num,DIM>>;

pub type CalculateBodyStateOutput<Num,const DIM:usize>=Sum<CalculatePositionStateOutput<Num,DIM>,CalculateAngularStateOutput<Num,DIM>>;

pub fn calculate_body_state<Num,const DIM:usize>(values : HToRef<CalculateBodyStateInput<Num,DIM>>)
->CalculateBodyStateOutput<Num,DIM>
where
	Num:RealField+Copy,
	Const<DIM>: DimNameToSoDimName + DimName,
	DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>+Allocator<DimNameToSoDimNameType<DIM>>,
    DimNameToSoDimNameType<DIM>:
        DimMin<DimNameToSoDimNameType<DIM>, Output = DimNameToSoDimNameType<DIM>>,

{
	let (a,b)=values.sculpt();
	calculate_position_state(a)+calculate_angular_state(b)
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
}

pub fn calculate_body_state_inplace<Num,const DIM:usize>(
	values:
	<HToRef<CalculateBodyStateInput<Num,DIM>> as Add<HToMut<CalculateBodyStateOutput<Num,DIM>>> >::Output)
where
	Num:RealField+Copy,
	Const<DIM>: DimNameToSoDimName + DimName,
	DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>+Allocator<DimNameToSoDimNameType<DIM>>,
    DimNameToSoDimNameType<DIM>:
        DimMin<DimNameToSoDimNameType<DIM>, Output = DimNameToSoDimNameType<DIM>>,
{
	let (inputs,outputs):(HToRef<CalculateBodyStateInput<Num,DIM>>,_)=values.sculpt();
	outputs.zip(calculate_body_state(inputs)).map(Poly(SetMut));
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
	calculate_body_state_inplace(hlist![mass,momentum,agi,agm,vel,kinetic,agv,agk]);
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