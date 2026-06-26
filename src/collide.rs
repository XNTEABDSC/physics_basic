
use nalgebra::RealField;

pub fn collide<Num:RealField+Copy>(a_mass:Num,a_vel:Num,b_mass:Num,b_vel:Num,e:Num)
->(Num,Num)
{
	let mass_sum=a_mass+b_mass;
	let v_dif=b_vel-a_vel;
	let a_mom=a_mass*a_vel;
	let b_mom=b_mass*b_vel;
	let a_v_2=(a_mom+b_mom + e*b_mass*(v_dif))/mass_sum;
	let b_v_2=(a_mom+b_mom + e*a_mass*(-v_dif))/mass_sum;
	(a_v_2,b_v_2)
}

pub fn collide_inelastic<Num:RealField+Copy>(a_mass:Num,a_vel:Num,b_mass:Num,b_vel:Num)
->Num
{
	(a_mass*a_vel+b_mass*b_vel)/(a_mass+b_mass)
}

