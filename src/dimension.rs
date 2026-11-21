use std::ops::{Add, Div, Mul, Rem, Sub};

use nalgebra::{
    Dim, DimMin, LU, OMatrix, OVector, RealField, SMatrix, SVector, ToConst, ToTypenum,
};
use num_traits::{Num, Zero};
use typenum::{U0, U1, U2};

pub use nalgebra::{DefaultAllocator,allocator::Allocator,Const,DimName};

/// N -> N*(N-1)/2
pub trait DimToSoDim:Dim
where
    Self: ToTypenum<
        Typenum: Sub<U1>
                     + Mul<
            <Self::Typenum as Sub<U1>>::Output,
            Output: Rem<U2, Output = U0> + Div<U2, Output: ToConst>,
        >,
    >,
{
    type SoDim:Dim+DimName;
    const SO_DIM: usize;
}

// type ToSoSizeTExp<N>= <Prod<N,Sum<N,N1>> as PartialDiv<P2>>::Output;

impl<const DIM: usize, TNum, SubTNum1, TNum2> DimToSoDim for Const<DIM>
where
    Self: ToTypenum<Typenum = TNum>,
    TNum: Sub<U1, Output = SubTNum1>
        + Mul<SubTNum1, Output: Div<U2, Output = TNum2> + Rem<U2, Output = U0>>,
    TNum2: ToConst,
{
    type SoDim = TNum2::Const;

    const SO_DIM: usize = TNum2::Const::DIM;
}

pub type DimToSoDimT<Dim> = <Dim as DimToSoDim>::SoDim;
pub type ConstDimToSoDimT<const DIM: usize> = <Const<DIM> as DimToSoDim>::SoDim;

pub trait DimSquare
where Self:DimMin<Self, Output = Self>{}

impl<T> DimSquare for T
where Self:DimMin<Self, Output = Self>{}

pub trait DimSquareSo
where 
	Self:DimToSoDim + DimName + DimSquare,{}

impl<const DIM:usize> DimSquareSo for Const<DIM>
where 
	Self:DimToSoDim + DimName + DimMin<Self, Output = Self>,{}


pub trait AllocatorSq<D:Dim> : Allocator<D,D>{}

impl<D:Dim,T> AllocatorSq<D> for T
where Self:Allocator<D,D>{}

pub trait AllocatorVM<D:Dim> : Allocator<D>+AllocatorSq<D>{}

impl<D:Dim,T> AllocatorVM<D> for T
where Self:Allocator<D>+AllocatorSq<D>{}

pub trait AllocatorDimVMSoDimVM<D>
where 
	D:DimSquareSo,
	<D as DimToSoDim>::SoDim:DimSquare,
	Self:AllocatorVM<D>+AllocatorVM<DimToSoDimT<D>>{}

impl<D:Dim,T> AllocatorDimVMSoDimVM<D> for T 
where 
	D:DimSquareSo,
	<D as DimToSoDim>::SoDim:DimSquare,
	Self:AllocatorVM<D>+AllocatorVM<DimToSoDimT<D>>{}


pub trait AllocatorSync<D:Dim,Num:RealField>:Allocator<D,Buffer<Num>:Send+Sync>{}
impl<D:Dim,Num:RealField,T> AllocatorSync<D,Num> for T
where Self:Allocator<D,Buffer<Num>:Send+Sync>{}

pub trait AllocatorSyncSq<D:Dim,Num:RealField> : Allocator<D,D,Buffer<Num>:Send+Sync>{}

impl<D:Dim,Num:RealField,T> AllocatorSyncSq<D,Num> for T
where Self:Allocator<D,D,Buffer<Num>:Send+Sync>{}

pub trait AllocatorSyncVMSq<D:Dim,Num:RealField> : AllocatorSync<D,Num>+AllocatorSyncSq<D,Num>{}

impl<D:Dim,Num:RealField,T> AllocatorSyncVMSq<D,Num> for T
where Self:AllocatorSync<D,Num>+AllocatorSyncSq<D,Num>{}

pub trait AllocatorSyncDimSqSoDimSq<D:DimToSoDim,Num:RealField>:AllocatorSyncVMSq<D,Num>+AllocatorSyncVMSq<DimToSoDimT<D>,Num>{}

impl<D:DimToSoDim,Num:RealField,T> AllocatorSyncDimSqSoDimSq<D,Num> for T
where T:AllocatorSyncVMSq<D,Num>+AllocatorSyncVMSq<DimToSoDimT<D>,Num>{}

/// N -> N*(N-1)/2
pub trait DimToDimAdd1:Dim
where
    Self: ToTypenum<
        Typenum: Add<U1, Output : ToConst>,
    >,
{
    type SoDim:Dim+DimName;
    const SO_DIM: usize;
}

// type ToSoSizeTExp<N>= <Prod<N,Sum<N,N1>> as PartialDiv<P2>>::Output;

impl<const DIM: usize, NumPlus1> DimToDimAdd1 for Const<DIM>
where
    Self: ToTypenum<
        Typenum: Add<U1, Output = NumPlus1>,
    >,
    NumPlus1:ToConst
{
    type SoDim = NumPlus1::Const;

    const SO_DIM: usize = NumPlus1::Const::DIM;
}


pub type DimToDimAdd1T<Dim> = <Dim as DimToDimAdd1>::SoDim;
pub type ConstDimToDimAdd1T<const DIM: usize> = <Const<DIM> as DimToDimAdd1>::SoDim;