use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

pub trait Ring:
    Sized
    + Clone
    + PartialEq
    + fmt::Display
    + Add<Output = Self>
    + AddAssign
    + for<'a> Add<&'a Self, Output = Self>
    + Sub<Output = Self>
    + SubAssign
    + for<'a> Sub<&'a Self, Output = Self>
    + Mul<Output = Self>
    + MulAssign
    + for<'a> Mul<&'a Self, Output = Self>
    + Neg<Output = Self>
    + One
    + Zero
{
    // A ring is an algebraic structure with addition and multiplication
    fn add_ref(&self, rhs: &Self) -> Self;
    fn sub_ref(&self, rhs: &Self) -> Self;
    fn mul_ref(&self, rhs: &Self) -> Self;

    fn add_assign_ref(&mut self, other: &Self);
    fn sub_assign_ref(&mut self, other: &Self);
    fn mul_assign_ref(&mut self, other: &Self);

    // Utility functions
    fn pow<M: Into<BigInt>>(&self, n: M) -> Self;
    fn get_value(&self) -> BigInt;
    fn from_value<M: Into<BigInt>>(value: M) -> Self;
    fn random_element(exclude_elements: &[Self]) -> Self;
}
