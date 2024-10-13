use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_traits::{One, Zero};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FieldElement {
    val: BigInt,
    k_modulus: BigInt,
}

impl FieldElement {
    // Constructor with optional k_modulus and generator_val.
    fn new(val: BigInt, k_modulus: Option<BigInt>) -> Self {
        let default_k_modulus = (3_i64 * (1 << 30) + 1).to_bigint().unwrap();
        let k_modulus = k_modulus.unwrap_or(default_k_modulus.clone());

        FieldElement {
            val: val % &k_modulus,
            k_modulus,
        }
    }

    // Zero element.
    pub fn zero(k_modulus: Option<BigInt>) -> Self {
        FieldElement::new(BigInt::zero(), k_modulus)
    }

    // Unit element.
    pub fn one(k_modulus: Option<BigInt>) -> Self {
        FieldElement::new(BigInt::one(), k_modulus)
    }

    // Typecasting from an integer or another FieldElement.
    fn typecast<T: Into<BigInt>>(other: T, k_modulus: &BigInt) -> FieldElement {
        let val: BigInt = other.into();
        FieldElement::new(val, Some(k_modulus.clone()))
    }

    // Inverse of the element using extended Euclidean algorithm.
    pub fn inverse(&self) -> Self {
        let (mut t, mut new_t) = (BigInt::zero(), BigInt::one());
        let (mut r, mut new_r) = (self.k_modulus.clone(), self.val.clone());
        while !new_r.is_zero() {
            let quotient = &r / &new_r;
            (t, new_t) = (new_t.clone(), &t - &quotient * &new_t);
            (r, new_r) = (new_r.clone(), &r - &quotient * &new_r);
        }
        assert_eq!(r, BigInt::one());
        FieldElement::new(t, Some(self.k_modulus.clone()))
    }

    // Check the order of an element.
    fn is_order(&self, n: u64) -> bool {
        assert!(n >= 1);
        let mut h = FieldElement::one(Some(self.k_modulus.clone()));
        for _ in 1..n {
            h = h * self.clone();
            if h == FieldElement::one(Some(self.k_modulus.clone())) {
                return false;
            }
        }
        h * self.clone() == FieldElement::one(Some(self.k_modulus.clone()))
    }

    // Serialize method.
    fn serialize(&self) -> String {
        self.val.to_string()
    }

    // Random element excluding a set of elements.
    fn random_element(k_modulus: BigInt, exclude_elements: &[FieldElement]) -> Self {
        let mut rng = rand::thread_rng();
        let mut fe = FieldElement::new(
            rng.gen_bigint_range(&BigInt::zero(), &k_modulus),
            Some(k_modulus.clone()),
        );
        while exclude_elements.contains(&fe) {
            fe = FieldElement::new(
                rng.gen_bigint_range(&BigInt::zero(), &k_modulus),
                Some(k_modulus.clone()),
            );
        }
        fe
    }
}

// Display trait implementation for pretty printing.
impl fmt::Display for FieldElement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let repr_val = (&self.val + &self.k_modulus / 2) % &self.k_modulus - &self.k_modulus / 2;
        write!(f, "{}", repr_val)
    }
}

// Arithmetic operations implementation for FieldElement.
impl Add for FieldElement {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        FieldElement::new(
            (&self.val + &other.val) % &self.k_modulus,
            Some(self.k_modulus.clone()),
        )
    }
}

impl Sub for FieldElement {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        FieldElement::new(
            (&self.val - &other.val) % &self.k_modulus,
            Some(self.k_modulus.clone()),
        )
    }
}

impl Mul for FieldElement {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        FieldElement::new(
            (&self.val * &other.val) % &self.k_modulus,
            Some(self.k_modulus.clone()),
        )
    }
}

impl Neg for FieldElement {
    type Output = Self;

    fn neg(self) -> Self {
        FieldElement::zero(Some(self.k_modulus.clone())) - self
    }
}

impl Div for FieldElement {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        self * other.inverse()
    }
}

// Power implementation for exponentiation.
impl FieldElement {
    fn pow(&self, mut n: u64) -> Self {
        assert!(n >= 0);
        let mut cur_pow = self.clone();
        let mut res = FieldElement::one(Some(self.k_modulus.clone()));
        while n > 0 {
            if n % 2 != 0 {
                res = res * cur_pow.clone();
            }
            n /= 2;
            cur_pow = cur_pow.clone() * cur_pow;
        }
        res
    }
}

// Implement conversion from BigInt for FieldElement.
impl From<BigInt> for FieldElement {
    fn from(val: BigInt) -> Self {
        FieldElement::new(val, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::ToBigInt;

    #[test]
    fn test_new_default_modulus() {
        let fe = FieldElement::new(10.to_bigint().unwrap(), None);
        assert_eq!(fe.val, 10.to_bigint().unwrap());
        assert_eq!(fe.k_modulus, (3_i64 * (1 << 30) + 1).to_bigint().unwrap());
    }

    #[test]
    fn test_new_custom_modulus() {
        let fe = FieldElement::new(10.to_bigint().unwrap(), Some(7.to_bigint().unwrap()));
        assert_eq!(fe.val, 10.to_bigint().unwrap() % 7.to_bigint().unwrap());
        assert_eq!(fe.k_modulus, 7.to_bigint().unwrap());
    }

    #[test]
    fn test_zero() {
        let fe_zero = FieldElement::zero(None);
        assert_eq!(fe_zero.val, BigInt::zero());
        assert_eq!(
            fe_zero.k_modulus,
            (3_i64 * (1 << 30) + 1).to_bigint().unwrap()
        );
    }

    #[test]
    fn test_one() {
        let fe_one = FieldElement::one(None);
        assert_eq!(fe_one.val, BigInt::one());
        assert_eq!(
            fe_one.k_modulus,
            (3_i64 * (1 << 30) + 1).to_bigint().unwrap()
        );
    }

    #[test]
    fn test_addition() {
        let fe1 = FieldElement::new(10.to_bigint().unwrap(), None);
        let fe2 = FieldElement::new(15.to_bigint().unwrap(), None);
        let fe3 = fe1 + fe2;
        assert_eq!(fe3.val, (10 + 15).to_bigint().unwrap() % fe3.k_modulus);
    }

    #[test]
    fn test_subtraction() {
        let fe1 = FieldElement::new(20.to_bigint().unwrap(), None);
        let fe2 = FieldElement::new(15.to_bigint().unwrap(), None);
        let fe3 = fe1 - fe2;
        assert_eq!(fe3.val, (20 - 15).to_bigint().unwrap() % fe3.k_modulus);
    }

    #[test]
    fn test_multiplication() {
        let fe1 = FieldElement::new(5.to_bigint().unwrap(), None);
        let fe2 = FieldElement::new(4.to_bigint().unwrap(), None);
        let fe3 = fe1 * fe2;
        assert_eq!(fe3.val, (5 * 4).to_bigint().unwrap() % fe3.k_modulus);
    }

    #[test]
    fn test_inverse() {
        let fe1 = FieldElement::new(7.to_bigint().unwrap(), Some(17.to_bigint().unwrap()));
        let fe_inv = fe1.inverse();
        assert_eq!(fe_inv.val, 5.to_bigint().unwrap());
        assert_eq!((fe1 * fe_inv).val, 1.to_bigint().unwrap());
    }

    #[test]
    fn test_negation() {
        let fe1 = FieldElement::new(10.to_bigint().unwrap(), None);
        let fe_neg = -fe1;
        assert_eq!(
            fe_neg.val,
            (0.to_bigint().unwrap() - 10.to_bigint().unwrap()) % fe_neg.k_modulus
        );
    }

    #[test]
    fn test_division() {
        let fe1 = FieldElement::new(12.to_bigint().unwrap(), Some(17.to_bigint().unwrap()));
        let fe2 = FieldElement::new(3.to_bigint().unwrap(), Some(17.to_bigint().unwrap()));
        let fe_div = fe1 / fe2.clone();
        assert_eq!(
            fe_div.val,
            (12.to_bigint().unwrap() * fe2.inverse().val) % fe_div.k_modulus
        );
    }

    #[test]
    fn test_exponentiation() {
        let fe1 = FieldElement::new(3.to_bigint().unwrap(), Some(17.to_bigint().unwrap()));
        let fe_exp = fe1.pow(4);
        assert_eq!(
            fe_exp.val,
            3.to_bigint().unwrap().pow(4) % 17.to_bigint().unwrap()
        );
    }

    #[test]
    fn test_is_order() {
        let fe1 = FieldElement::new(2.to_bigint().unwrap(), Some(17.to_bigint().unwrap()));
        assert!(fe1.is_order(8));
        assert!(!fe1.is_order(3));
    }

    #[test]
    fn test_random_element() {
        let excluded_elements = vec![
            FieldElement::new(2.to_bigint().unwrap(), None),
            FieldElement::new(3.to_bigint().unwrap(), None),
        ];
        let fe_random = FieldElement::random_element(17.to_bigint().unwrap(), &excluded_elements);
        assert!(!excluded_elements.contains(&fe_random));
    }
}
