use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_traits::{One, Zero};
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
    + One
    + Zero
{
    // A ring is an algebraic structure with addition and multiplication

    // Utility functions
    fn pow(&self, n: BigInt) -> Self;
    fn get_value(&self) -> BigInt;
    fn from_value<M: Into<BigInt>>(value: M) -> Self;
    fn random_element(exclude_elements: &[Self]) -> Self;
}
