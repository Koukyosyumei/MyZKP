use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_traits::{One, Signed, Zero};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

const DEFAULT_K_MODULES: i64 = 3_i64 * (1 << 30) + 1;

pub trait Field:
    Sized
    + Clone
    + PartialEq
    + fmt::Display
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
{
    // A commutative ring with a multiplicative identity element
    // where every non-zero element has a multiplicative inverse is called a field.
    fn one(modulus: Option<BigInt>) -> Self;
    fn inverse(&self) -> Self;

    // Utility functions
    fn zero(modulus: Option<BigInt>) -> Self;
    fn mul_scalar(&self, n: i64) -> Self;
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FiniteFieldElement {
    value: BigInt,
    modulus: BigInt,
}

impl FiniteFieldElement {
    // Constructor with optional modulus and generator_value.
    fn new(value: BigInt, modulus: Option<BigInt>) -> Self {
        let default_modulus = (DEFAULT_K_MODULES).to_bigint().unwrap();
        let modulus = modulus.unwrap_or(default_modulus);

        let mut value_sanitized = value % &modulus;
        if value_sanitized < 0_i32.to_bigint().unwrap() {
            value_sanitized += modulus.clone();
        }

        FiniteFieldElement {
            value: value_sanitized,
            modulus,
        }
    }

    // Check the order of an element.
    fn is_order(&self, n: u64) -> bool {
        assert!(n >= 1);
        let identity = FiniteFieldElement::one(Some(self.modulus.clone()));
        let mut h = identity.clone();
        for _ in 1..n {
            h = h * self.clone();
            if h == identity {
                return false;
            }
        }
        h * self.clone() == identity
    }

    // Serialize method.
    fn serialize(&self) -> String {
        self.value.to_string()
    }

    // Random element excluding a set of elements.
    pub fn random_element(exclude_elements: &[FiniteFieldElement]) -> Self {
        let modulus = (DEFAULT_K_MODULES).to_bigint().unwrap();
        let mut rng = rand::thread_rng();
        let mut fe = FiniteFieldElement::new(
            rng.gen_bigint_range(&BigInt::zero(), &modulus),
            Some(modulus.clone()),
        );
        while exclude_elements.contains(&fe) {
            fe = FiniteFieldElement::new(
                rng.gen_bigint_range(&BigInt::zero(), &modulus),
                Some(modulus.clone()),
            );
        }
        fe
    }
}

impl Field for FiniteFieldElement {
    fn zero(modulus: Option<BigInt>) -> Self {
        FiniteFieldElement::new(BigInt::zero(), modulus)
    }

    fn one(modulus: Option<BigInt>) -> Self {
        FiniteFieldElement::new(BigInt::one(), modulus)
    }

    fn mul_scalar(&self, n: i64) -> Self {
        FiniteFieldElement::new(
            (&self.value * n.to_bigint().unwrap()) % &self.modulus,
            Some(self.modulus.clone()),
        )
    }

    fn inverse(&self) -> Self {
        let mut t = BigInt::zero();
        let mut new_t = BigInt::one();
        let mut r = self.modulus.clone();
        let mut new_r = self.value.clone();

        // Ensure r and new_r are non-negative
        if r.is_negative() {
            r += &self.modulus;
        }
        if new_r.is_negative() {
            new_r += &self.modulus;
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
        t %= &self.modulus;
        if t.is_negative() {
            t += &self.modulus;
        }

        FiniteFieldElement::new(t, Some(self.modulus.clone()))
    }
}

// Display trait implementation for pretty printing.
impl fmt::Display for FiniteFieldElement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let repr_value = (&self.value + &self.modulus / 2) % &self.modulus - &self.modulus / 2;
        write!(f, "{}", repr_value)
    }
}

// Arithmetic operations implementation for FiniteFieldElement.
impl Add for FiniteFieldElement {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        FiniteFieldElement::new(
            (&self.value + &other.value) % &self.modulus,
            Some(self.modulus.clone()),
        )
    }
}

impl Sub for FiniteFieldElement {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        FiniteFieldElement::new(
            (&self.value - &other.value) % &self.modulus,
            Some(self.modulus.clone()),
        )
    }
}

impl Mul for FiniteFieldElement {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        FiniteFieldElement::new(
            (&self.value * &other.value) % &self.modulus,
            Some(self.modulus.clone()),
        )
    }
}

impl Neg for FiniteFieldElement {
    type Output = Self;

    fn neg(self) -> Self {
        FiniteFieldElement::zero(Some(self.modulus.clone())) - self
    }
}

impl Div for FiniteFieldElement {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        self * other.inverse()
    }
}

impl FiniteFieldElement {
    pub fn pow_scalar(&self, mut n: u64) -> Self {
        assert!(n >= 0);
        let mut cur_pow = self.clone();
        let mut res = FiniteFieldElement::one(Some(self.modulus.clone()));
        while n > 0 {
            if n % 2 != 0 {
                res = res * cur_pow.clone();
            }
            n /= 2;
            cur_pow = cur_pow.clone() * cur_pow;
        }
        res
    }

    pub fn pow(&self, n: FiniteFieldElement) -> Self {
        assert!(n.value.is_positive());
        let mut cur_pow = self.clone();
        let mut res = FiniteFieldElement::one(Some(self.modulus.clone()));
        let mut n_val = n.value;
        while n_val.is_positive() && !n_val.is_zero() {
            if n_val.clone() % 2 != BigInt::zero() {
                res = res * cur_pow.clone();
            }
            n_val = n_val / 2_i64.to_bigint().unwrap();
            cur_pow = cur_pow.clone() * cur_pow;
        }
        res
    }
}

// Implement conversion from BigInt for FiniteFieldElement.
impl From<BigInt> for FiniteFieldElement {
    fn from(value: BigInt) -> Self {
        FiniteFieldElement::new(value, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::ToBigInt;

    #[test]
    fn test_new_default_modulus() {
        let fe = FiniteFieldElement::new(10.to_bigint().unwrap(), None);
        assert_eq!(fe.value, 10.to_bigint().unwrap());
        assert_eq!(fe.modulus, (3_i64 * (1 << 30) + 1).to_bigint().unwrap());
    }

    #[test]
    fn test_new_custom_modulus() {
        let fe = FiniteFieldElement::new(10.to_bigint().unwrap(), Some(7.to_bigint().unwrap()));
        assert_eq!(fe.value, 10.to_bigint().unwrap() % 7.to_bigint().unwrap());
        assert_eq!(fe.modulus, 7.to_bigint().unwrap());
    }

    #[test]
    fn test_zero() {
        let fe_zero = FiniteFieldElement::zero(None);
        assert_eq!(fe_zero.value, BigInt::zero());
        assert_eq!(
            fe_zero.modulus,
            (3_i64 * (1 << 30) + 1).to_bigint().unwrap()
        );
    }

    #[test]
    fn test_one() {
        let fe_one = FiniteFieldElement::one(None);
        assert_eq!(fe_one.value, BigInt::one());
        assert_eq!(fe_one.modulus, (3_i64 * (1 << 30) + 1).to_bigint().unwrap());
    }

    #[test]
    fn test_addition() {
        let fe1 = FiniteFieldElement::new(10.to_bigint().unwrap(), None);
        let fe2 = FiniteFieldElement::new(15.to_bigint().unwrap(), None);
        let fe3 = fe1 + fe2;
        assert_eq!(fe3.value, (10 + 15).to_bigint().unwrap() % fe3.modulus);
    }

    #[test]
    fn test_subtraction() {
        let fe1 = FiniteFieldElement::new(20.to_bigint().unwrap(), None);
        let fe2 = FiniteFieldElement::new(15.to_bigint().unwrap(), None);
        let fe3 = fe1 - fe2;
        assert_eq!(fe3.value, (20 - 15).to_bigint().unwrap() % fe3.modulus);
    }

    #[test]
    fn test_multiplication() {
        let fe1 = FiniteFieldElement::new(5.to_bigint().unwrap(), None);
        let fe2 = FiniteFieldElement::new(4.to_bigint().unwrap(), None);
        let fe3 = fe1 * fe2;
        assert_eq!(fe3.value, (5 * 4).to_bigint().unwrap() % fe3.modulus);
    }

    #[test]
    fn test_inverse() {
        let fe1 = FiniteFieldElement::new(7.to_bigint().unwrap(), Some(17.to_bigint().unwrap()));
        let fe_inv = fe1.inverse();
        assert_eq!(fe_inv.value, 5.to_bigint().unwrap());
        assert_eq!((fe1 * fe_inv).value, 1.to_bigint().unwrap());
    }

    #[test]
    fn test_negation() {
        let fe1 = FiniteFieldElement::new(10.to_bigint().unwrap(), None);
        let fe_neg = -fe1;
        assert_eq!(
            fe_neg.value,
            ((-10 + DEFAULT_K_MODULES).to_bigint().unwrap()) % fe_neg.modulus
        );
    }

    #[test]
    fn test_division() {
        let fe1 = FiniteFieldElement::new(12.to_bigint().unwrap(), Some(17.to_bigint().unwrap()));
        let fe2 = FiniteFieldElement::new(3.to_bigint().unwrap(), Some(17.to_bigint().unwrap()));
        let fe_div = fe1 / fe2.clone();
        assert_eq!(
            fe_div.value,
            (12.to_bigint().unwrap() * fe2.inverse().value) % fe_div.modulus
        );
    }

    #[test]
    fn test_exponentiation() {
        let fe1 = FiniteFieldElement::new(3.to_bigint().unwrap(), Some(17.to_bigint().unwrap()));
        let fe_exp = fe1.pow_scalar(4);
        assert_eq!(
            fe_exp.value,
            3.to_bigint().unwrap().pow(4) % 17.to_bigint().unwrap()
        );
    }

    #[test]
    fn test_is_order() {
        let fe1 = FiniteFieldElement::new(2.to_bigint().unwrap(), Some(17.to_bigint().unwrap()));
        assert!(fe1.is_order(8));
        assert!(!fe1.is_order(3));
    }

    #[test]
    fn test_random_element() {
        let excluded_elements = vec![
            FiniteFieldElement::new(2.to_bigint().unwrap(), None),
            FiniteFieldElement::new(3.to_bigint().unwrap(), None),
        ];
        let fe_random = FiniteFieldElement::random_element(&excluded_elements);
        assert!(!excluded_elements.contains(&fe_random));
    }
}
