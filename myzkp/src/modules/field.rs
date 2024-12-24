//! # Finite Field Element Implementation
//!
//! This module provides a generic implementation of finite field elements and related operations.
//! It includes traits and structures for working with finite fields, as well as utility functions
//! for arithmetic operations in finite fields.
//!
//! ## Key Components
//!
//! - `Field` trait: Defines the interface for field operations.
//! - `FiniteFieldElement<M>` struct: Represents an element in a finite field.
//! - `ModulusValue` trait: Defines the modulus for a specific finite field.
//!
//! ## Features
//!
//! - Arithmetic operations: addition, subtraction, multiplication, division, and exponentiation.
//! - Modular arithmetic with arbitrary precision integers.
//! - Random element generation.
//! - Order checking for field elements.
//! - Macro for easy definition of new modulus types.
//!
//! ## Usage
//!
//! To use this module, you need to define a modulus type using the `define_myzkp_modulus_type` macro,
//! and then create `FiniteFieldElement` instances with that modulus.
//!
//! ```
//! use std::str::FromStr;
//! use paste::paste;
//! use num_bigint::BigInt;
//! use lazy_static::lazy_static;
//! use myzkp::define_myzkp_modulus_type;
//! use myzkp::modules::ring::Ring;
//! use myzkp::modules::field::ModulusValue;
//! use myzkp::modules::field::FiniteFieldElement;
//!
//! define_myzkp_modulus_type!(Mod17, "17");
//! let a = FiniteFieldElement::<Mod17>::from_value(5);
//! let b = FiniteFieldElement::<Mod17>::from_value(3);
//! let c = a + b; // Performs addition modulo 17
//! ```
//!
//! ## Note
//!
//! This implementation uses the `num_bigint` crate for arbitrary-precision integer arithmetic,
//! allowing for operations on large finite fields commonly used in cryptographic applications.

use std::fmt;
use std::fmt::Debug;
use std::hash::Hash;
use std::hash::Hasher;
use std::marker::PhantomData;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::str::FromStr;

use lazy_static::lazy_static;
use num_bigint::{BigInt, RandBigInt};
use num_traits::{One, Signed, Zero};
use paste::paste;

use crate::modules::ring::Ring;
use crate::modules::utils::extended_euclidean;

/// Trait representing a field in abstract algebra.
///
/// This trait extends the `Ring` trait and adds operations specific to fields,
/// such as division and finding multiplicative inverses.
pub trait Field: Ring + Div<Output = Self> + PartialEq + Eq + Hash {
    /// Computes the multiplicative inverse of the element.
    fn inverse(&self) -> Self;
    fn div_ref(&self, other: &Self) -> Self;
}

/// Represents an element in a finite field.
///
/// The type parameter `M` specifies the modulus of the field.
#[derive(Debug, Clone)]
pub struct FiniteFieldElement<M> {
    pub value: BigInt,
    _phantom: PhantomData<M>,
}

/// Trait for defining the modulus of a finite field.
pub trait ModulusValue: Debug + Clone + Hash {
    fn modulus() -> &'static BigInt;
}

impl<M: ModulusValue> FiniteFieldElement<M> {
    /// Creates a new `FiniteFieldElement` with the given value.
    ///
    /// The value is automatically reduced modulo the field's modulus.
    pub fn new(value: BigInt) -> Self {
        let modulus = M::modulus();
        let value_sanitized = value % modulus;

        FiniteFieldElement {
            value: value_sanitized,
            _phantom: PhantomData,
        }
    }

    /// Ensures the element's value is within the correct range for the field.
    pub fn sanitize(&self) -> Self {
        let modulus = M::modulus();
        let mut value_sanitized = &self.value % modulus;
        if value_sanitized < BigInt::zero() {
            value_sanitized += modulus;
        }
        FiniteFieldElement {
            value: value_sanitized,
            _phantom: PhantomData,
        }
    }

    /// Checks if the element has the given order in the field.
    ///
    /// # Arguments
    ///
    /// * `n` - The order to check.
    ///
    /// # Returns
    ///
    /// `true` if the element has order `n`, `false` otherwise.
    pub fn is_order(&self, n: u64) -> bool {
        assert!(n >= 1);
        let identity = FiniteFieldElement::one();
        let mut h = identity.clone();
        for _ in 1..n {
            h = h * self;
            if h == identity {
                return false;
            }
        }
        h * self == identity
    }
}

impl<M: ModulusValue> Zero for FiniteFieldElement<M> {
    fn zero() -> Self {
        FiniteFieldElement::<M>::new(BigInt::zero())
    }

    fn is_zero(&self) -> bool {
        self.value.is_zero()
    }
}

impl<M: ModulusValue> One for FiniteFieldElement<M> {
    fn one() -> Self {
        FiniteFieldElement::<M>::new(BigInt::one())
    }
}

impl<M: ModulusValue> Hash for FiniteFieldElement<M> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.value.hash(state);
    }
}

impl<M: ModulusValue> Field for FiniteFieldElement<M> {
    fn inverse(&self) -> Self {
        let modulus = M::modulus();
        let mut r = modulus.clone();
        let mut new_r = self.value.clone();

        // Ensure r and new_r are non-negative
        if r.is_negative() {
            r += modulus;
        }
        if new_r.is_negative() {
            new_r += modulus;
        }

        let (_, _, mut t) = extended_euclidean(r, new_r);

        // At this point, r should be 1 if the inverse exists
        //if final_r == BigInt::one() {
        // Ensure t is within the correct range
        t %= modulus;
        if t.is_negative() {
            t += modulus;
        }

        FiniteFieldElement::<M>::new(t)
        //} else {
        //    panic!("r={}: Inverse does not exist", r)
        //}
    }

    fn div_ref(&self, other: &Self) -> Self {
        self.mul_ref(&other.inverse())
    }
}

impl<M: ModulusValue> Ring for FiniteFieldElement<M> {
    fn add_ref(&self, other: &Self) -> Self {
        FiniteFieldElement::<M>::new(&self.value + &other.value)
    }

    fn mul_ref(&self, other: &Self) -> Self {
        FiniteFieldElement::<M>::new(&self.value * &other.value)
    }

    fn sub_ref(&self, other: &Self) -> Self {
        FiniteFieldElement::<M>::new(&self.value - &other.value)
    }

    fn pow<V: Into<BigInt>>(&self, n: V) -> Self {
        let mut exponent: BigInt = n.into();
        let mut base = self.clone();

        if exponent.is_negative() {
            // x^{-n} = (x^{-1})^{n}
            base = base.inverse();
            exponent = -exponent;
        }

        let mut result = FiniteFieldElement::one();
        while exponent.is_positive() && !exponent.is_zero() {
            if &exponent % 2 != BigInt::zero() {
                result = result * &base;
            }
            exponent /= 2;
            base = base.mul_ref(&base);
        }
        result
    }

    fn from_value<V: Into<BigInt>>(value: V) -> Self {
        FiniteFieldElement::<M>::new(value.into())
    }

    fn get_value(&self) -> BigInt {
        self.value.clone()
    }

    // Random element excluding a set of elements.
    fn random_element(exclude_elements: &[FiniteFieldElement<M>]) -> Self {
        let modulus = M::modulus();
        let mut rng = rand::thread_rng();
        let mut fe = FiniteFieldElement::<M>::new(rng.gen_bigint_range(&BigInt::zero(), modulus));
        while exclude_elements.contains(&fe) {
            fe = FiniteFieldElement::<M>::new(rng.gen_bigint_range(&BigInt::zero(), modulus));
        }
        fe
    }
}

// Display trait implementation for pretty printing.
impl<M: ModulusValue> fmt::Display for FiniteFieldElement<M> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        //let modulus = M::modulus();
        //let repr_value = (&self.value + modulus / 2) % modulus - modulus / 2;
        write!(f, "{}", &self.value)
    }
}

impl<M: ModulusValue> PartialEq for FiniteFieldElement<M> {
    fn eq(&self, other: &Self) -> bool {
        self.sanitize().value == other.sanitize().value
    }
}

impl<M: ModulusValue> Eq for FiniteFieldElement<M> {}

// Arithmetic operations implementation for FiniteFieldElement<M>.
impl<M: ModulusValue> Add for FiniteFieldElement<M> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        self.add_ref(&other)
    }
}

impl<'a, M: ModulusValue> Add<&'a FiniteFieldElement<M>> for FiniteFieldElement<M> {
    type Output = FiniteFieldElement<M>;

    fn add(self, other: &'a FiniteFieldElement<M>) -> FiniteFieldElement<M> {
        self.add_ref(other)
    }
}

impl<M: ModulusValue> Sub for FiniteFieldElement<M> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self.sub_ref(&other)
    }
}

impl<'a, M: ModulusValue> Sub<&'a FiniteFieldElement<M>> for FiniteFieldElement<M> {
    type Output = Self;

    fn sub(self, other: &'a FiniteFieldElement<M>) -> FiniteFieldElement<M> {
        self.sub_ref(&other)
    }
}

impl<M: ModulusValue> Mul<FiniteFieldElement<M>> for FiniteFieldElement<M> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        self.mul_ref(&other)
    }
}

impl<'a, M: ModulusValue> Mul<&'a FiniteFieldElement<M>> for FiniteFieldElement<M> {
    type Output = FiniteFieldElement<M>;

    fn mul(self, other: &'a FiniteFieldElement<M>) -> FiniteFieldElement<M> {
        self.mul_ref(other)
    }
}

impl<M: ModulusValue> Neg for FiniteFieldElement<M> {
    type Output = Self;

    fn neg(self) -> Self {
        FiniteFieldElement::<M>::new(BigInt::zero()) - self
    }
}

impl<M: ModulusValue> Div for FiniteFieldElement<M> {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        self.div_ref(&other)
    }
}

/// Macro for defining a new modulus type.
///
/// This macro creates a new type that implements `ModulusValue` with the specified modulus.
///
/// # Example
///
/// ```
/// use std::str::FromStr;
/// use paste::paste;
/// use num_bigint::BigInt;
/// use lazy_static::lazy_static;
/// use myzkp::define_myzkp_modulus_type;
/// use myzkp::modules::ring::Ring;
/// use myzkp::modules::field::ModulusValue;
/// use myzkp::modules::field::FiniteFieldElement;
///
/// define_myzkp_modulus_type!(Mod17, "17");
/// ```
#[macro_export]
macro_rules! define_myzkp_modulus_type {
    ($name:ident, $modulus:expr) => {
        paste! {#[derive(Debug, Hash, Clone)]
            pub struct $name;

            lazy_static! {
                static ref [<MODULUS_ $name>]: BigInt = BigInt::from_str($modulus).unwrap();
            }

            impl ModulusValue for $name {
                fn modulus() -> &'static BigInt {
                    &[<MODULUS_ $name>]
                }
            }
        }
    };
}

// Pre-defined modulus for EIP-197
define_myzkp_modulus_type!(
    ModEIP197,
    "21888242871839275222246405745257275088548364400416034343698204186575808495617"
);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modules::ring::Ring;
    use num_bigint::ToBigInt;

    define_myzkp_modulus_type!(Mod7, "7");
    define_myzkp_modulus_type!(Mod17, "17");
    define_myzkp_modulus_type!(Mod31, "31");

    #[test]
    fn test_new_default_modulus() {
        let fe = FiniteFieldElement::<ModEIP197>::from_value(10);
        assert_eq!(fe.value, 10.to_bigint().unwrap());
    }

    #[test]
    fn test_new_custom_modulus() {
        let fe = FiniteFieldElement::<Mod7>::from_value(10);
        assert_eq!(fe.value, 10.to_bigint().unwrap() % 7.to_bigint().unwrap());
    }

    #[test]
    fn test_zero() {
        let fe_zero = FiniteFieldElement::<ModEIP197>::zero();
        assert_eq!(fe_zero.value, BigInt::zero());
    }

    #[test]
    fn test_one() {
        let fe_one = FiniteFieldElement::<ModEIP197>::one();
        assert_eq!(fe_one.value, BigInt::one());
    }

    #[test]
    fn test_addition() {
        let fe1 = FiniteFieldElement::<ModEIP197>::from_value(10);
        let fe2 = FiniteFieldElement::<ModEIP197>::from_value(15);
        let fe3 = fe1 + fe2;
        assert_eq!(fe3.value, (10 + 15).to_bigint().unwrap());
    }

    #[test]
    fn test_subtraction() {
        let fe1 = FiniteFieldElement::<ModEIP197>::from_value(20);
        let fe2 = FiniteFieldElement::<ModEIP197>::from_value(15);
        let fe3 = fe1 - fe2;
        assert_eq!(fe3.value, (20 - 15).to_bigint().unwrap());
    }

    #[test]
    fn test_multiplication() {
        let fe1 = FiniteFieldElement::<ModEIP197>::from_value(5);
        let fe2 = FiniteFieldElement::<ModEIP197>::from_value(4);
        let fe3 = fe1 * fe2;
        assert_eq!(fe3.value, (5 * 4).to_bigint().unwrap());
    }

    #[test]
    fn test_inverse() {
        let fe1 = FiniteFieldElement::<Mod17>::from_value(7);
        let fe_inv = fe1.inverse();
        assert_eq!(fe_inv.value, 5.to_bigint().unwrap());
        assert_eq!((fe1 * fe_inv).value, 1.to_bigint().unwrap());
    }

    #[test]
    fn test_negation() {
        let fe1 = FiniteFieldElement::<ModEIP197>::from_value(10);
        let fe_neg = -fe1;
        assert_eq!(fe_neg.value, ((-10_i128).to_bigint().unwrap()));
    }

    #[test]
    fn test_division() {
        let fe1 = FiniteFieldElement::<Mod17>::from_value(12);
        let fe2 = FiniteFieldElement::<Mod17>::from_value(3);
        let fe_div = fe1.div_ref(&fe2);
        assert_eq!(
            fe_div.value,
            (12.to_bigint().unwrap() * fe2.inverse().value) % 17
        );
    }

    #[test]
    fn test_exponentiation() {
        let fe1 = FiniteFieldElement::<Mod17>::from_value(3);
        let fe_exp = fe1.pow(4);
        assert_eq!(
            fe_exp.value,
            3.to_bigint().unwrap().pow(4) % 17.to_bigint().unwrap()
        );
    }

    #[test]
    fn test_is_order() {
        let fe1 = FiniteFieldElement::<Mod17>::from_value(2);
        assert!(fe1.is_order(8));
        assert!(!fe1.is_order(3));
    }

    #[test]
    fn test_random_element() {
        let excluded_elements = vec![
            FiniteFieldElement::<ModEIP197>::from_value(2),
            FiniteFieldElement::<ModEIP197>::from_value(3),
        ];
        let fe_random = FiniteFieldElement::random_element(&excluded_elements);
        assert!(!excluded_elements.contains(&fe_random));
    }

    #[test]
    fn test_eq() {
        assert_eq!(
            FiniteFieldElement::<Mod31>::from_value(-23),
            FiniteFieldElement::<Mod31>::from_value(8)
        )
    }
}
