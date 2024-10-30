use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_traits::{One, Signed, Zero};
use std::fmt;
use std::fmt::Debug;
use std::hash::Hash;
use std::hash::Hasher;
use std::marker::PhantomData;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::str::FromStr;

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

#[derive(Debug, Clone)]
pub struct FiniteFieldElement<M> {
    pub value: BigInt,
    _phantom: PhantomData<M>,
}

pub trait ModulusValue: Debug + Clone + Hash {
    fn modulus() -> BigInt;
}

impl<M: ModulusValue> FiniteFieldElement<M> {
    // Constructor with optional modulus and generator_value.
    pub fn new(value: BigInt) -> Self {
        let modulus = M::modulus();
        let value_sanitized = value % &modulus;

        //if value_sanitized < 0_i32.to_bigint().unwrap() {
        //    value_sanitized += modulus.clone();
        //}

        FiniteFieldElement {
            value: value_sanitized,
            _phantom: PhantomData,
        }
    }

    pub fn sanitize(&self) -> Self {
        let modulus = M::modulus();
        let mut value_sanitized = self.value.clone() % &modulus;
        if value_sanitized < 0_i32.to_bigint().unwrap() {
            value_sanitized += modulus.clone();
        }
        FiniteFieldElement {
            value: value_sanitized,
            _phantom: PhantomData,
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
impl<M: ModulusValue> Ring for FiniteFieldElement<M> {
    fn zero() -> Self {
        FiniteFieldElement::<M>::new(BigInt::zero())
    }

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
    fn get_value(&self) -> BigInt {
        self.value.clone()
    }
    fn inverse(&self) -> Self {
        let modulus = M::modulus();
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

        FiniteFieldElement::<M>::new(t)
    }

    fn pow(&self, n: BigInt) -> Self {
        let mut exponent = n;
        let mut base = self.clone();

        if exponent.is_negative() {
            // x^{-n} = (x^{-1})^{n}
            base = base.inverse();
            exponent = -exponent;
        }

        let mut result = FiniteFieldElement::one();
        while exponent.is_positive() && !exponent.is_zero() {
            if exponent.clone() % 2 != BigInt::zero() {
                result = result * base.clone();
            }
            exponent /= 2;
            base = base.clone() * base;
        }
        result
    }

    fn from_value<V: Into<BigInt>>(value: V) -> Self {
        FiniteFieldElement::<M>::new(value.into())
    }

    // Random element excluding a set of elements.
    fn random_element(exclude_elements: &[FiniteFieldElement<M>]) -> Self {
        let modulus = M::modulus();
        let mut rng = rand::thread_rng();
        let mut fe = FiniteFieldElement::<M>::new(rng.gen_bigint_range(&BigInt::zero(), &modulus));
        while exclude_elements.contains(&fe) {
            fe = FiniteFieldElement::<M>::new(rng.gen_bigint_range(&BigInt::zero(), &modulus));
        }
        fe
    }
}

// Display trait implementation for pretty printing.
impl<M: ModulusValue> fmt::Display for FiniteFieldElement<M> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        //let modulus = M::modulus();
        //let repr_value = (&self.value + &modulus / 2) % &modulus - &modulus / 2;
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
        let modulus = M::modulus();
        FiniteFieldElement::<M>::new((&self.value + &other.value) % &modulus)
    }
}

impl<M: ModulusValue> Sub for FiniteFieldElement<M> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let modulus = M::modulus();
        FiniteFieldElement::<M>::new((&self.value - &other.value) % &modulus)
    }
}

impl<M: ModulusValue> Mul<FiniteFieldElement<M>> for FiniteFieldElement<M> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let modulus = M::modulus();
        FiniteFieldElement::<M>::new((&self.value * &other.value) % &modulus)
    }
}

impl<M: ModulusValue> Mul<i64> for FiniteFieldElement<M> {
    type Output = Self;

    fn mul(self, n: i64) -> Self {
        let modulus = M::modulus();
        FiniteFieldElement::<M>::new((&self.value * n.to_bigint().unwrap()) % &modulus)
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
        self * other.inverse()
    }
}

#[derive(Debug, Hash, Clone)]
pub struct ModDEFAULT;
impl ModulusValue for ModDEFAULT {
    fn modulus() -> BigInt {
        BigInt::from_str(
            "21888242871839275222246405745257275088548364400416034343698204186575808495617",
        )
        .unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::ToBigInt;

    #[test]
    fn test_new_default_modulus() {
        let fe = FiniteFieldElement::<ModDEFAULT>::from_value(10);
        assert_eq!(fe.value, 10.to_bigint().unwrap());
    }

    #[derive(Debug, Hash, Clone)]
    struct Mod7;
    impl ModulusValue for Mod7 {
        fn modulus() -> BigInt {
            BigInt::from(7)
        }
    }

    #[test]
    fn test_new_custom_modulus() {
        let fe = FiniteFieldElement::<Mod7>::from_value(10);
        assert_eq!(fe.value, 10.to_bigint().unwrap() % 7.to_bigint().unwrap());
    }

    #[test]
    fn test_zero() {
        let fe_zero = FiniteFieldElement::<ModDEFAULT>::zero();
        assert_eq!(fe_zero.value, BigInt::zero());
    }

    #[test]
    fn test_one() {
        let fe_one = FiniteFieldElement::<ModDEFAULT>::one();
        assert_eq!(fe_one.value, BigInt::one());
    }

    #[test]
    fn test_addition() {
        let fe1 = FiniteFieldElement::<ModDEFAULT>::from_value(10);
        let fe2 = FiniteFieldElement::<ModDEFAULT>::from_value(15);
        let fe3 = fe1 + fe2;
        assert_eq!(fe3.value, (10 + 15).to_bigint().unwrap());
    }

    #[test]
    fn test_subtraction() {
        let fe1 = FiniteFieldElement::<ModDEFAULT>::from_value(20);
        let fe2 = FiniteFieldElement::<ModDEFAULT>::from_value(15);
        let fe3 = fe1 - fe2;
        assert_eq!(fe3.value, (20 - 15).to_bigint().unwrap());
    }

    #[test]
    fn test_multiplication() {
        let fe1 = FiniteFieldElement::<ModDEFAULT>::from_value(5);
        let fe2 = FiniteFieldElement::<ModDEFAULT>::from_value(4);
        let fe3 = fe1 * fe2;
        assert_eq!(fe3.value, (5 * 4).to_bigint().unwrap());
    }

    #[derive(Debug, Hash, Clone)]
    struct Mod17;
    impl ModulusValue for Mod17 {
        fn modulus() -> BigInt {
            BigInt::from(17)
        }
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
        let fe1 = FiniteFieldElement::<ModDEFAULT>::from_value(10);
        let fe_neg = -fe1;
        assert_eq!(
            fe_neg.value,
            ((-10_i128).to_bigint().unwrap()) % DEFAULT_K_MODULES
        );
    }

    #[test]
    fn test_division() {
        let fe1 = FiniteFieldElement::<Mod17>::from_value(12);
        let fe2 = FiniteFieldElement::<Mod17>::from_value(3);
        let fe_div = fe1 / fe2.clone();
        assert_eq!(
            fe_div.value,
            (12.to_bigint().unwrap() * fe2.inverse().value) % 17
        );
    }

    #[test]
    fn test_exponentiation() {
        let fe1 = FiniteFieldElement::<Mod17>::from_value(3);
        let fe_exp = fe1.pow(4.to_bigint().unwrap());
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
            FiniteFieldElement::<ModDEFAULT>::from_value(2),
            FiniteFieldElement::<ModDEFAULT>::from_value(3),
        ];
        let fe_random = FiniteFieldElement::random_element(&excluded_elements);
        assert!(!excluded_elements.contains(&fe_random));
    }

    #[derive(Debug, Hash, Clone)]
    struct Mod31;
    impl ModulusValue for Mod31 {
        fn modulus() -> BigInt {
            BigInt::from(31)
        }
    }

    #[test]
    fn test_eq() {
        assert_eq!(
            FiniteFieldElement::<Mod31>::from_value(-23),
            FiniteFieldElement::<Mod31>::from_value(8)
        )
    }
}
