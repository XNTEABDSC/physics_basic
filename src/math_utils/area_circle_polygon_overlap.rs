use std::cell::LazyCell;

//pub type Num=fixed::types::I32F32;
use crate::num::*;
use cordic::{atan, atan2, sqrt};

#[inline]
pub fn area_circle_triangle<'a>(r:Num,r_sq:Num,d:Num,y1:Num,y1_angle:&'a LazyCell<Num,impl FnOnce() -> Num>,y2:Num,y2_angle:&'a LazyCell<Num,impl FnOnce() -> Num>)->Num{
    if d.is_zero(){
        return Num::ZERO;
    }
    if d.is_negative(){
        return area_circle_triangle_inner(r,r_sq,d.wrapping_neg(),
        y1.wrapping_neg(),
        &LazyCell::new( ||{(**y1_angle) -Num::FRAC_PI_2} ),
        y2.wrapping_neg(),
        &LazyCell::new( ||{(**y2_angle) -Num::FRAC_PI_2} ));
    }
    return area_circle_triangle_inner(r,r_sq,d,y1,y1_angle,y2,y2_angle);
}

/// Computes the overlap area between a circle centered at (0,0) and a triangle
/// with vertices (0,0), (d, y1), and (d, y2).
/// - `r`: Radius of the circle.
/// - `r_sq`: Precomputed radius squared (r * r).
/// - `r_sq_over2`: Precomputed r_sq / 2 for efficiency.
/// - `d`: X-coordinate of the vertical leg of the triangle.
/// - `y1`, `y2`: Y-coordinates defining the triangle's vertical segment at x = d.
/// Returns the signed overlap area; negative if y1 < y2 and the orientation is reversed.
fn area_circle_triangle_inner<'a>(r:Num,r_sq:Num,d:Num,y1:Num,y1_angle:&'a LazyCell<Num,impl FnOnce() -> Num>,y2:Num,y2_angle:&'a LazyCell<Num,impl FnOnce() -> Num>)->Num{
    let  r_sq_over2=r_sq.wrapping_shr(1);
    if d>=r{
        return r_sq_over2*(**y2_angle-**y1_angle);
    }
    
    let d_sq=d*d;
    let y_d=sqrt( r_sq-d_sq );
    let area_d=Num::ONE.wrapping_shr(1)*d*y_d;
    let atan_frac_y_d_d=LazyCell::new(||atan(y_d/d));

    struct YData<'a,Awawdawdaw:FnOnce() -> Num>{
        pub y_sign_neg:bool,
        pub y_abs_greater_than_y_d:bool,
        pub y_abs:Num,
        pub y_angle:&'a LazyCell<Num,Awawdawdaw>
    }
    let new_ydata1=|y:Num,y_angle:&'a LazyCell<Num,_>|{
        let y_sign_neg=y.is_negative();
        let y_abs=y.abs();
        let y_abs_greater_than_y_d=y_abs>y_d;
        YData{
            //y,
            y_sign_neg,
            y_abs,
            y_angle,
            y_abs_greater_than_y_d,
        }
    };
    let new_ydata2=|y:Num,y_angle:&'a LazyCell<Num,_>|{
        let y_sign_neg=y.is_negative();
        let y_abs=y.abs();
        let y_abs_greater_than_y_d=y_abs>y_d;
        YData{
            //y,
            y_sign_neg,
            y_abs,
            y_angle,
            y_abs_greater_than_y_d,
        }
    };
    let y1data=new_ydata1(y1,y1_angle);
    let y2data=new_ydata2(y2,y2_angle);

    
    if y1data.y_abs_greater_than_y_d && y2data.y_abs_greater_than_y_d && (y1data.y_sign_neg==y2data.y_sign_neg){
        return r_sq_over2*(**y2_angle-**y1_angle);
    } 

    let helper_f1=|ydata:YData<_>|{

        let mut area:Num;
        if ydata.y_abs_greater_than_y_d{
            area=r_sq_over2*((**ydata.y_angle).abs()-*atan_frac_y_d_d) + area_d;
        }else {
            area=Num::ONE.wrapping_shr(1)*d*ydata.y_abs;
        }
        if ydata.y_sign_neg{
            area=area.wrapping_neg();
        }
        return area;
    };
    let helper_f2=|ydata:YData<_>|{

        let mut area:Num;
        if ydata.y_abs_greater_than_y_d{
            area=r_sq_over2*((**ydata.y_angle).abs()-*atan_frac_y_d_d) + area_d;
        }else {
            area=Num::ONE.wrapping_shr(1)*d*ydata.y_abs;
        }
        if ydata.y_sign_neg{
            area=area.wrapping_neg();
        }
        return area;
    };

    let y1_area=helper_f1(y1data);
    let y2_area=helper_f2(y2data);

    return  y2_area.wrapping_sub(y1_area);
}

#[test]
fn test_area_circle_triangle() {
    fn test_area_circle_triangle_fn(r:Num,r_sq:Num,d:Num,y1:Num,y2:Num)
    {
        let y1_angle=LazyCell::new(||atan2(y1, d));
        let y2_angle=LazyCell::new(||atan2(y2, d));
        let a=area_circle_triangle(r, r_sq, d, y1, &y1_angle, y2, &y2_angle);
        let b=area_circle_triangle_old(r, r_sq, d, y1, y2);
        assert_eq!(a,b,"test {r}, {r_sq}, {d}, {y1}, {y2}")
    }

    test_area_circle_triangle_fn(2.into(),4.into(),1.into(),(-1).into(),1.into());
    test_area_circle_triangle_fn(2.into(),4.into(),2.into(),(-1).into(),1.into());
    test_area_circle_triangle_fn(2.into(),4.into(),1.into(),(-3).into(),3.into());
    test_area_circle_triangle_fn(2.into(),4.into(),1.into(),(3).into(),(-3).into());
    test_area_circle_triangle_fn(2.into(),4.into(),1.into(),(-1).into(),3.into());
    test_area_circle_triangle_fn(2.into(),4.into(),1.into(),(3).into(),4.into());

}


pub fn area_circle_rectangle(r:Num,r_sq:Num,x1:Num,y1:Num,x2:Num,y2:Num)->Num {
    let atan_1=LazyCell::new(||{atan2(y1,x1)});
    let atan_2=LazyCell::new(||{atan2(y1,x2)});
    let atan_3=LazyCell::new(||{atan2(y2,x2)});
    let atan_4=LazyCell::new(||{atan2(y2,x1)});
    fn atan_lc_rotatel90(v:&LazyCell<Num,impl FnOnce()->Num>)->LazyCell<Num,impl FnOnce()->Num>{
        return LazyCell::new(||{**v+Num::FRAC_PI_2});
    }
    fn atan_lc_rotater90(v:&LazyCell<Num,impl FnOnce()->Num>)->LazyCell<Num,impl FnOnce()->Num>{
        return LazyCell::new(||{**v-Num::FRAC_PI_2});
    }
    let atan_1_rl90=atan_lc_rotatel90(&atan_1);
    let atan_2_rl90=atan_lc_rotatel90(&atan_2);

    let atan_3_rr90=atan_lc_rotater90(&atan_3);
    let atan_4_rr90=atan_lc_rotater90(&atan_4);
    

    let area1=
    area_circle_triangle(r, r_sq, y1.wrapping_neg(), x1,&atan_1_rl90, x2,&atan_2_rl90);
    let area2=
    area_circle_triangle(r, r_sq, x2, y1,&atan_2, y2,&atan_3);
    let area3=
    area_circle_triangle(r, r_sq, y2, x2.wrapping_neg(),&atan_3_rr90, x1.wrapping_neg(),&atan_4_rr90);
    let area4=
    //area_circle_triangle_old(r, r_sq, x1.wrapping_neg(), y2.wrapping_neg(), y1.wrapping_neg());
    area_circle_triangle(r, r_sq, x1, y2,&atan_4, y1,&atan_1);
    return area1+area2+area3+area4;
}

#[test]
fn test_area_circle_rectangle() {
    let r = Num::ONE * 2; // r = 2
    let r_sq = r * r;

    // Case 1: Rectangle fully inside circle
    assert_eq!(
        area_circle_rectangle(r, r_sq, Num::ZERO, Num::ZERO, Num::ONE, Num::ONE),
        Num::ONE,
        "Rectangle [0,0] to [1,1], r=2, expected area 1"
    );

    // Case 2: Rectangle fully outside circle
    let r_small = Num::ONE; // r = 1
    let r_sq_small = r_small * r_small;
    assert_eq!(
        area_circle_rectangle(
            r_small,
            r_sq_small,
            Num::ONE * 2,
            Num::ONE * 2,
            Num::ONE * 3,
            Num::ONE * 3
        ),
        Num::ZERO,
        "Rectangle [2,2] to [3,3], r=1, expected area 0"
    );

    // Case 3: Rectangle over origin, fully inside circle
    assert_eq!(
        area_circle_rectangle(
            r,
            r_sq,
            Num::ONE.wrapping_neg(),
            Num::ONE.wrapping_neg(),
            Num::ONE,
            Num::ONE
        ),
        Num::ONE * 4,
        "Rectangle [-1,-1] to [1,1], r=2, expected area 4"
    );

    // Case 4: Circle fully inside rectangle
    let r_small = Num::ONE; // r = 1
    let r_sq_small = r_small * r_small;
    assert_eq!(
        area_circle_rectangle(
            r_small,
            r_sq_small,
            Num::ONE * -2,
            Num::ONE * -2,
            Num::ONE * 2,
            Num::ONE * 2
        ),
        Num::PI * r_small * r_small,
        "Rectangle [-2,-2] to [2,2], r=1, expected area π"
    );

    // Case 5: Partial overlap
    let partial_area = area_circle_rectangle(r, r_sq, Num::ONE, Num::ONE, Num::ONE * 2, Num::ONE * 2);
    assert!(
        partial_area > Num::ZERO,
        "Rectangle [1,1] to [2,2], r=2, expected area > 0"
    );

    // Case 6: Rectangle touching circle boundary
    assert_eq!(
        area_circle_rectangle(r, r_sq, Num::ONE * 2, Num::ZERO, Num::ONE * 3, Num::ONE),
        Num::ZERO,
        "Rectangle [2,0] to [3,1], r=2, expected area 0"
    );

    // Case 7: Rectangle intersecting circle in one quadrant
    let overlap_area = area_circle_rectangle(r, r_sq, Num::ONE, Num::ONE, Num::ONE * 3, Num::ONE * 3);
    assert!(
        overlap_area > Num::ZERO && overlap_area < Num::ONE * 4,
        "Rectangle [1,1] to [3,3], r=2, expected area between 0 and 4"
    );
}

/// Computes the overlap area between a circle centered at (0,0) and a triangle
/// with vertices (0,0), (d, y1), and (d, y2).
/// - `r`: Radius of the circle.
/// - `r_sq`: Precomputed radius squared (r * r).
/// - `r_sq_over2`: Precomputed r_sq / 2 for efficiency.
/// - `d`: X-coordinate of the vertical leg of the triangle.
/// - `y1`, `y2`: Y-coordinates defining the triangle's vertical segment at x = d.
/// Returns the signed overlap area; negative if y1 < y2 and the orientation is reversed.
fn area_circle_triangle_old(r:Num,r_sq:Num,mut d:Num,mut y1:Num,mut y2:Num)->Num{
    
    if d.is_zero(){
        return Num::ZERO;
    }
    
    if d.is_negative(){
        d=d.wrapping_neg();
        y2=y2.wrapping_neg();
        y1=y1.wrapping_neg();
    }else {

    }
    let  r_sq_over2=r_sq.wrapping_shr(1);
    if d>=r{
        return r_sq_over2*(atan(y2/d)-atan(y1/d));
    }
    else {
        let d_sq=d*d;
        let y_d=sqrt( r_sq-d_sq );
        let area_d=Num::ONE.wrapping_shr(1)*d*y_d;
        let atan_frac_y_d_d=atan(y_d/d);

        let y1_sign_positive;
        if y1.is_negative(){
            y1_sign_positive=false;
            y1=y1.wrapping_neg();
        }else{
            y1_sign_positive=true;
        }

        let y2_sign_neg;
        if y2.is_negative(){
            y2_sign_neg=true;
            y2=y2.wrapping_neg();
        }else {
            y2_sign_neg=false;
        }

        let y1_greater_than_y_d=y1>y_d;
        let y2_greater_than_y_d=y2>y_d;

        if y1_greater_than_y_d&&y2_greater_than_y_d && (y1_sign_positive!=y2_sign_neg){
            let mut y2_angle=atan(y2/d);
            let mut y1_angle=atan(y1/d);
            if y1_sign_positive{
                y1_angle=y1_angle.wrapping_neg();
            }
            if y2_sign_neg{
                y2_angle=y2_angle.wrapping_neg();
            }
            return r_sq_over2*(y2_angle+y1_angle);
        }

        let helper_f=|y:Num,y_greater_than_y_d,sign_neg:bool|{
            let mut area;
            if y_greater_than_y_d{
                area=r_sq_over2*(atan(y2/d)-atan_frac_y_d_d) + area_d;
            }else {
                area=Num::ONE.wrapping_shr(1)*d*y;
            }
            if sign_neg{
                area=area.wrapping_neg();
            }
            return area;
        };

        let y2_area=helper_f(y2,y2_greater_than_y_d,y2_sign_neg);
        let y1_area_neg=helper_f(y1,y1_greater_than_y_d,y1_sign_positive);

        return  y2_area.wrapping_add(y1_area_neg) ;

    }

}

#[test]
fn test_area_circle_triangle_old(){
    let r=Num::ONE*2;
    let r_sq=r*r;
    assert_eq!(
        area_circle_triangle_old(r,r_sq,1.into(),(-1).into(),1.into()),
        Num::ONE
    );
    assert_eq!(
        area_circle_triangle_old(r,r_sq,2.into(),(-2).into(),2.into()),
        Num::ONE*2*Num::PI/2
    );
    assert_eq!(
        area_circle_triangle_old(r,r_sq,2.into(),2.into(),(-2).into()),
        -Num::ONE*2*Num::PI/2
    );
    assert_eq!(
        area_circle_triangle_old(r,r_sq,1.into(),(-1).into(),2.into()),
        Num::ONE/2
        +Num::ONE/2*Num::sqrt(3.into())
        +r*(
            cordic::atan::<Num>(2.into())
            -atan((Num::ONE*3).sqrt())
        )
    );

}
/*
    let r_sq=r*r;
    let d_sq=d*d;
    let y_d;
    let sign_neg=false;
    if y<0{
        sign=true;
        y=-y;
    }
    if d<r {
        y_d=sqrt(r_sq-d_sq)
    }else{
        y_d=0
    }
    let area;
    if y>y_d{
        area=r*( atan(y/d)-atan(y_d/d) )+0.5*d*y_d;
    }
    else{
        area=0.5*d*y;
    }
    if y>0{
    }
    return area
 */

fn area_circle_rectangle_old(r:Num,r_sq:Num,x1:Num,y1:Num,x2:Num,y2:Num)->Num {
    let area1=area_circle_triangle_old(r, r_sq, y1.wrapping_neg(), x1, x2);
    let area2=area_circle_triangle_old(r, r_sq, x2, y1, y2);
    let area3=area_circle_triangle_old(r, r_sq, y2, x2.wrapping_neg(), x1.wrapping_neg());
    let area4=area_circle_triangle_old(r, r_sq, x1.wrapping_neg(), y2.wrapping_neg(), y1.wrapping_neg());
    return area1+area2+area3+area4;
}

#[test]
fn test_area_circle_rectangle_old() {
    let r = Num::ONE * 2; // r = 2
    let r_sq = r * r;

    // Case 1: Rectangle fully inside circle
    assert_eq!(
        area_circle_rectangle_old(r, r_sq, Num::ZERO, Num::ZERO, Num::ONE, Num::ONE),
        Num::ONE,
        "Rectangle [0,0] to [1,1], r=2, expected area 1"
    );

    // Case 2: Rectangle fully outside circle
    let r_small = Num::ONE; // r = 1
    let r_sq_small = r_small * r_small;
    assert_eq!(
        area_circle_rectangle_old(
            r_small,
            r_sq_small,
            Num::ONE * 2,
            Num::ONE * 2,
            Num::ONE * 3,
            Num::ONE * 3
        ),
        Num::ZERO,
        "Rectangle [2,2] to [3,3], r=1, expected area 0"
    );

    // Case 3: Rectangle over origin, fully inside circle
    assert_eq!(
        area_circle_rectangle_old(
            r,
            r_sq,
            Num::ONE.wrapping_neg(),
            Num::ONE.wrapping_neg(),
            Num::ONE,
            Num::ONE
        ),
        Num::ONE * 4,
        "Rectangle [-1,-1] to [1,1], r=2, expected area 4"
    );

    // Case 4: Circle fully inside rectangle
    let r_small = Num::ONE; // r = 1
    let r_sq_small = r_small * r_small;
    assert_eq!(
        area_circle_rectangle_old(
            r_small,
            r_sq_small,
            Num::ONE * -2,
            Num::ONE * -2,
            Num::ONE * 2,
            Num::ONE * 2
        ),
        Num::PI * r_small * r_small,
        "Rectangle [-2,-2] to [2,2], r=1, expected area π"
    );

    // Case 5: Partial overlap
    assert!(
        area_circle_rectangle_old(r, r_sq, Num::ONE, Num::ONE, Num::ONE * 2, Num::ONE * 2) > Num::ZERO,
        "Rectangle [1,1] to [2,2], r=2, expected area > 0"
    );
}