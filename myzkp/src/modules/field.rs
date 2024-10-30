use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_traits::{One, Signed, Zero};
use std::fmt;
use std::hash::Hash;
use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::modules::ring::Ring;

pub const DEFAULT_K_MODULES: i128 = 3_i128 * (1 << 30) + 1;

pub trait Field: Ring + Div<Output = Self> + PartialEq + Eq + Hash {
    // A commutative ring with a multiplicative identity element
    // where every non-zero element has a multiplicative inverse is called a field.
    fn inverse(&self) -> Self;

    // Utility functions
    fn get_value(&self) -> BigInt;
    fn pow(&self, n: BigInt) -> Self;
    fn from_value<M: Into<BigInt>>(value: M) -> Self;
    fn random_element(exclude_elements: &[Self]) -> Self;
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct FiniteFieldElement<const MODULUS: i128> {
    pub value: BigInt,
}

impl<const MODULUS: i128> FiniteFieldElement<MODULUS> {
    // Constructor with optional modulus and generator_value.
    pub fn new(value: BigInt) -> Self {
        let modulus = MODULUS.to_bigint().unwrap();
        let mut value_sanitized = value % &modulus;
        if value_sanitized < 0_i32.to_bigint().unwrap() {
            value_sanitized += modulus.clone();
        }

        FiniteFieldElement {
            value: value_sanitized,
        }
    }

    // Check the order of an element.
    pub fn is_order(&self, n: u64) -> bool {
        assert!(n >= 1);
        let identity = FiniteFieldElement::one();
        let mut h = identity.clone();
        for _ in 1..n {
            h = h * self.clone();
            if h == identity {
                return false;
            }
        }
        h * self.clone() == identity
    }

    /*
    // Serialize method.
    fn serialize(&self) -> String {
        self.value.to_string()
    }
    */
}
impl<const MODULUS: i128> Ring for FiniteFieldElement<MODULUS> {
    fn zero() -> Self {
        FiniteFieldElement::<MODULUS>::new(BigInt::zero())
    }

    fn one() -> Self {
        FiniteFieldElement::<MODULUS>::new(BigInt::one())
    }
}

impl<const MODULUS: i128> Field for FiniteFieldElement<MODULUS> {
    fn get_value(&self) -> BigInt {
        self.value.clone()
    }
    fn inverse(&self) -> Self {
        let modulus = MODULUS.to_bigint().unwrap();
        let mut t = BigInt::zero();
        let mut new_t = BigInt::one();
        let mut r = modulus.clone();
        let mut new_r = self.value.clone();

        // Ensure r and new_r are non-negative
        if r.is_negative() {
            r += &modulus;
        }
        if new_r.is_negative() {
            new_r += &modulus;
        }

        while !new_r.is_zero() {
            let quotient = &r / &new_r;

            // Update (t, new_t) = (new_t, t - quotient * new_t)
            let temp_t = new_t.clone();
            new_t = &t - &quotient * &new_t;
            t = temp_t;

            // Update (r, new_r) = (new_r, r - quotient * new_r)
            let temp_r = new_r.clone();
            new_r = &r - &quotient * &new_r;
            r = temp_r;
        }

        // At this point, r should be 1 if the inverse exists
        assert_eq!(r, BigInt::one(), "Inverse does not exist");

        // Ensure t is within the correct range
        t %= &modulus;
        if t.is_negative() {
            t += &modulus;
        }

        FiniteFieldElement::<MODULUS>::new(t)
    }

    fn pow(&self, mut n: BigInt) -> Self {
        assert!(!n.is_negative());
        let mut cur_pow = self.clone();
        let mut res = FiniteFieldElement::one();
        while n.is_positive() && !n.is_zero() {
            if n.clone() % 2 != BigInt::zero() {
                res = res * cur_pow.clone();
            }
            n /= 2;
            cur_pow = cur_pow.clone() * cur_pow;
        }
        res
    }

    fn from_value<M: Into<BigInt>>(value: M) -> Self {
        FiniteFieldElement::<MODULUS>::new(value.into())
    }

    // Random element excluding a set of elements.
    fn random_element(exclude_elements: &[FiniteFieldElement<MODULUS>]) -> Self {
        let modulus = MODULUS.to_bigint().unwrap();
        let mut rng = rand::thread_rng();
        let mut fe =
            FiniteFieldElement::<MODULUS>::new(rng.gen_bigint_range(&BigInt::zero(), &modulus));
        while exclude_elements.contains(&fe) {
            fe =
                FiniteFieldElement::<MODULUS>::new(rng.gen_bigint_range(&BigInt::zero(), &modulus));
        }
        fe
    }
}

// Display trait implementation for pretty printing.
impl<const MODULUS: i128> fmt::Display for FiniteFieldElement<MODULUS> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let modulus = MODULUS.to_bigint().unwrap();
        let repr_value = (&self.value + &modulus / 2) % &modulus - &modulus / 2;
        write!(f, "{}", repr_value)
    }
}

// Arithmetic operations implementation for FiniteFieldElement<MODULUS>.
impl<const MODULUS: i128> Add for FiniteFieldElement<MODULUS> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let modulus = MODULUS.to_bigint().unwrap();
        FiniteFieldElement::<MODULUS>::new((&self.value + &other.value) % &modulus)
    }
}

impl<const MODULUS: i128> Sub for FiniteFieldElement<MODULUS> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let modulus = MODULUS.to_bigint().unwrap();
        FiniteFieldElement::<MODULUS>::new((&self.value - &other.value) % &modulus)
    }
}

impl<const MODULUS: i128> Mul<FiniteFieldElement<MODULUS>> for FiniteFieldElement<MODULUS> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let modulus = MODULUS.to_bigint().unwrap();
        FiniteFieldElement::<MODULUS>::new((&self.value * &other.value) % &modulus)
    }
}

impl<const MODULUS: i128> Mul<i64> for FiniteFieldElement<MODULUS> {
    type Output = Self;

    fn mul(self, n: i64) -> Self {
        let modulus = MODULUS.to_bigint().unwrap();
        FiniteFieldElement::<MODULUS>::new((&self.value * n.to_bigint().unwrap()) % &modulus)
    }
}

impl<const MODULUS: i128> Neg for FiniteFieldElement<MODULUS> {
    type Output = Self;

    fn neg(self) -> Self {
        FiniteFieldElement::<MODULUS>::new(BigInt::zero()) - self
    }
}

impl<const MODULUS: i128> Div for FiniteFieldElement<MODULUS> {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        self * other.inverse()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::ToBigInt;

    #[test]
    fn test_new_default_modulus() {
        let fe = FiniteFieldElement::<DEFAULT_K_MODULES>::new(10.to_bigint().unwrap());
        assert_eq!(fe.value, 10.to_bigint().unwrap());
    }

    #[test]
    fn test_new_custom_modulus() {
        let fe = FiniteFieldElement::<7_i128>::new(10.to_bigint().unwrap());
        assert_eq!(fe.value, 10.to_bigint().unwrap() % 7.to_bigint().unwrap());
    }

    #[test]
    fn test_zero() {
        let fe_zero = FiniteFieldElement::<DEFAULT_K_MODULES>::zero();
        assert_eq!(fe_zero.value, BigInt::zero());
    }

    #[test]
    fn test_one() {
        let fe_one = FiniteFieldElement::<DEFAULT_K_MODULES>::one();
        assert_eq!(fe_one.value, BigInt::one());
    }

    #[test]
    fn test_addition() {
        let fe1 = FiniteFieldElement::<DEFAULT_K_MODULES>::new(10.to_bigint().unwrap());
        let fe2 = FiniteFieldElement::<DEFAULT_K_MODULES>::new(15.to_bigint().unwrap());
        let fe3 = fe1 + fe2;
        assert_eq!(
            fe3.value,
            (10 + 15).to_bigint().unwrap() % DEFAULT_K_MODULES
        );
    }

    #[test]
    fn test_subtraction() {
        let fe1 = FiniteFieldElement::<DEFAULT_K_MODULES>::new(20.to_bigint().unwrap());
        let fe2 = FiniteFieldElement::<DEFAULT_K_MODULES>::new(15.to_bigint().unwrap());
        let fe3 = fe1 - fe2;
        assert_eq!(
            fe3.value,
            (20 - 15).to_bigint().unwrap() % DEFAULT_K_MODULES
        );
    }

    #[test]
    fn test_multiplication() {
        let fe1 = FiniteFieldElement::<DEFAULT_K_MODULES>::new(5.to_bigint().unwrap());
        let fe2 = FiniteFieldElement::<DEFAULT_K_MODULES>::new(4.to_bigint().unwrap());
        let fe3 = fe1 * fe2;
        assert_eq!(fe3.value, (5 * 4).to_bigint().unwrap() % DEFAULT_K_MODULES);
    }

    #[test]
    fn test_inverse() {
        let fe1 = FiniteFieldElement::<17>::new(7.to_bigint().unwrap());
        let fe_inv = fe1.inverse();
        assert_eq!(fe_inv.value, 5.to_bigint().unwrap());
        assert_eq!((fe1 * fe_inv).value, 1.to_bigint().unwrap());
    }

    #[test]
    fn test_negation() {
        let fe1 = FiniteFieldElement::<DEFAULT_K_MODULES>::new(10.to_bigint().unwrap());
        let fe_neg = -fe1;
        assert_eq!(
            fe_neg.value,
            ((-10_i128 + DEFAULT_K_MODULES).to_bigint().unwrap()) % DEFAULT_K_MODULES
        );
    }

    #[test]
    fn test_division() {
        let fe1 = FiniteFieldElement::<17>::new(12.to_bigint().unwrap());
        let fe2 = FiniteFieldElement::<17>::new(3.to_bigint().unwrap());
        let fe_div = fe1 / fe2.clone();
        assert_eq!(
            fe_div.value,
            (12.to_bigint().unwrap() * fe2.inverse().value) % 17
        );
    }

    #[test]
    fn test_exponentiation() {
        let fe1 = FiniteFieldElement::<17>::new(3.to_bigint().unwrap());
        let fe_exp = fe1.pow(4.to_bigint().unwrap());
        assert_eq!(
            fe_exp.value,
            3.to_bigint().unwrap().pow(4) % 17.to_bigint().unwrap()
        );
    }

    #[test]
    fn test_is_order() {
        let fe1 = FiniteFieldElement::<17>::new(2.to_bigint().unwrap());
        assert!(fe1.is_order(8));
        assert!(!fe1.is_order(3));
    }

    #[test]
    fn test_random_element() {
        let excluded_elements = vec![
            FiniteFieldElement::<DEFAULT_K_MODULES>::new(2.to_bigint().unwrap()),
            FiniteFieldElement::<DEFAULT_K_MODULES>::new(3.to_bigint().unwrap()),
        ];
        let fe_random = FiniteFieldElement::random_element(&excluded_elements);
        assert!(!excluded_elements.contains(&fe_random));
    }
}
