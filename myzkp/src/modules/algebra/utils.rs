use std::ops::{Add, BitAnd, Div, Mul, Rem, Shl, ShrAssign, Sub};

use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

/// Computes the greatest common divisor (GCD) of two numbers `a` and `b`
/// using the Extended Euclidean Algorithm, along with the coefficients
/// `s` and `t` that satisfy the equation:
///
/// `GCD(a, b) = s * a + t * b`
///
/// This function works with any type `F` that supports arithmetic operations
/// and implements the necessary traits.
///
/// # Arguments
///
/// - `a`: The first number.
/// - `b`: The second number.
///
/// # Returns
///
/// A tuple `(gcd, s, t)` where:
/// - `gcd`: The greatest common divisor of `a` and `b`.
/// - `s`: The coefficient such that `s * a + t * b = gcd`.
/// - `t`: The coefficient such that `s * a + t * b = gcd`.
///
/// # Type Constraints
///
/// The type `F` must implement the following traits:
/// - `Clone`, `PartialEq`, `Zero`, `One`
/// - Arithmetic operations: `Sub`, `Div`, `Rem`
/// - Reference-based arithmetic operations: `&F - &F`, `&F / &F`, `&F % &F`, `&F * &F`
///
/// # Examples
///
/// ```
/// use num_bigint::BigInt;
/// use num_traits::{One, Zero};
/// use myzkp::modules::algebra::utils::extended_euclidean;
///
/// fn main() {
///     let a = BigInt::from(56);
///     let b = BigInt::from(15);
///
///     // Compute the GCD and coefficients
///     let (gcd, s, t) = extended_euclidean(a.clone(), b.clone());
///
///     // Verify the result
///     assert_eq!(gcd, BigInt::from(1)); // GCD of 56 and 15 is 1
///     assert_eq!(gcd, &s * &a + &t * &b); // s * a + t * b = GCD
/// }
/// ```
pub fn extended_euclidean<F>(a: F, b: F) -> (F, F, F)
where
    F: Clone + PartialEq + Zero + One,
    for<'a> &'a F: Sub<&'a F, Output = F>
        + Div<&'a F, Output = F>
        + Rem<&'a F, Output = F>
        + Mul<&'a F, Output = F>,
{
    let mut r0 = a;
    let mut r1 = b;
    let mut s0 = F::one();
    let mut s1 = F::zero();
    let mut t0 = F::zero();
    let mut t1 = F::one();

    while !r1.is_zero() {
        let q = &r0 / &r1;
        let r = &r0 % &r1;
        r0 = r1;
        r1 = r;
        let new_s = &s0 - &(&q * &s1);
        s0 = s1;
        s1 = new_s;
        let new_t = &t0 - &(&q * &t1);
        t0 = t1;
        t1 = new_t;
    }

    (r0, s0, t0)
}

pub fn mod_inverse<F>(a: &F, modulus: &F) -> F
where
    F: Clone + PartialEq + Zero + One + Signed,
    for<'a> &'a F: Add<&'a F, Output = F>
        + Sub<&'a F, Output = F>
        + Div<&'a F, Output = F>
        + Rem<&'a F, Output = F>
        + Mul<&'a F, Output = F>,
{
    let r = modulus.clone();
    let mut new_r = a.clone();

    if new_r.is_negative() {
        new_r = &new_r + modulus;
    }

    let (_, _, mut t) = extended_euclidean(r, new_r);
    t = &t % modulus;
    if t.is_negative() {
        t = &t + modulus;
    }

    t
}

pub fn mod_pow<F>(x: &F, y: &F, modulus: &F) -> F
where
    F: Clone + PartialEq + Zero + One + Signed + Shl<u32, Output = F> + ShrAssign<i32>,
    for<'a> &'a F: Add<&'a F, Output = F>
        + Sub<&'a F, Output = F>
        + Div<&'a F, Output = F>
        + Rem<&'a F, Output = F>
        + Mul<&'a F, Output = F>
        + BitAnd<&'a F, Output = F>,
{
    let mut base = x.clone();
    let mut exponent = y.clone();

    if exponent.is_negative() {
        base = mod_inverse(&base, modulus);
        exponent = -exponent;
    }
    let one = F::one();

    base = &base % modulus;
    let mut result = F::one();
    while exponent.is_positive() && !exponent.is_zero() {
        if (&exponent & &one) == one {
            result = &(&result * &base) % modulus;
        }
        exponent >>= 1;
        base = &(&base * &base) % modulus;
    }
    result
}
