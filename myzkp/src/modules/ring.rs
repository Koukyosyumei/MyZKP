use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

pub trait Ring:
    Sized
    + Clone
    + PartialEq
    + fmt::Display
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
{
    // A ring is an algebraic structure with addition and multiplication
    fn zero() -> Self;
    fn one() -> Self;
}
