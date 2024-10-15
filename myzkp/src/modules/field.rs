use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_traits::{One, Signed, Zero};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ModuloFieldElement {
    value: BigInt,
    modulus: BigInt,
}

const DEFAULT_K_MODULES: i64 = 3_i64 * (1 << 30) + 1;

impl ModuloFieldElement {
    // Constructor with optional modulus and generator_value.
    fn new(value: BigInt, modulus: Option<BigInt>) -> Self {
        let default_modulus = (DEFAULT_K_MODULES).to_bigint().unwrap();
        let modulus = modulus.unwrap_or(default_modulus);

        let mut value_sanitized = value % &modulus;
        if value_sanitized < 0_i32.to_bigint().unwrap() {
            value_sanitized += modulus.clone();
        }

        ModuloFieldElement {
            value: value_sanitized,
            modulus,
        }
    }

    // Zero element.
    pub fn zero(modulus: Option<BigInt>) -> Self {
        ModuloFieldElement::new(BigInt::zero(), modulus)
    }

    // Unit element.
    pub fn one(modulus: Option<BigInt>) -> Self {
        ModuloFieldElement::new(BigInt::one(), modulus)
    }

    // Typecasting from an integer or another ModuloFieldElement.
    fn typecast<T: Into<BigInt>>(other: T, modulus: &BigInt) -> ModuloFieldElement {
        let value: BigInt = other.into();
        ModuloFieldElement::new(value, Some(modulus.clone()))
    }

    // Inverse of the element using extended Euclidean algorithm.
    pub fn inverse(&self) -> Self {
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

        ModuloFieldElement::new(t, Some(self.modulus.clone()))
    }

    // Check the order of an element.
    fn is_order(&self, n: u64) -> bool {
        assert!(n >= 1);
        let identity = ModuloFieldElement::one(Some(self.modulus.clone()));
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
    pub fn random_element(exclude_elements: &[ModuloFieldElement]) -> Self {
        let modulus = (DEFAULT_K_MODULES).to_bigint().unwrap();
        let mut rng = rand::thread_rng();
        let mut fe = ModuloFieldElement::new(
            rng.gen_bigint_range(&BigInt::zero(), &modulus),
            Some(modulus.clone()),
        );
        while exclude_elements.contains(&fe) {
            fe = ModuloFieldElement::new(
                rng.gen_bigint_range(&BigInt::zero(), &modulus),
                Some(modulus.clone()),
            );
        }
        fe
    }
}

// Display trait implementation for pretty printing.
impl fmt::Display for ModuloFieldElement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let repr_value = (&self.value + &self.modulus / 2) % &self.modulus - &self.modulus / 2;
        write!(f, "{}", repr_value)
    }
}

// Arithmetic operations implementation for ModuloFieldElement.
impl Add for ModuloFieldElement {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        ModuloFieldElement::new(
            (&self.value + &other.value) % &self.modulus,
            Some(self.modulus.clone()),
        )
    }
}

impl Sub for ModuloFieldElement {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        ModuloFieldElement::new(
            (&self.value - &other.value) % &self.modulus,
            Some(self.modulus.clone()),
        )
    }
}

impl Mul for ModuloFieldElement {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        ModuloFieldElement::new(
            (&self.value * &other.value) % &self.modulus,
            Some(self.modulus.clone()),
        )
    }
}

impl Neg for ModuloFieldElement {
    type Output = Self;

    fn neg(self) -> Self {
        ModuloFieldElement::zero(Some(self.modulus.clone())) - self
    }
}

impl Div for ModuloFieldElement {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        self * other.inverse()
    }
}

impl ModuloFieldElement {
    pub fn mul_scalar(&self, n: i64) -> Self {
        ModuloFieldElement::new(
            (&self.value * n.to_bigint().unwrap()) % &self.modulus,
            Some(self.modulus.clone()),
        )
    }

    pub fn pow_scalar(&self, mut n: u64) -> Self {
        assert!(n >= 0);
        let mut cur_pow = self.clone();
        let mut res = ModuloFieldElement::one(Some(self.modulus.clone()));
        while n > 0 {
            if n % 2 != 0 {
                res = res * cur_pow.clone();
            }
            n /= 2;
            cur_pow = cur_pow.clone() * cur_pow;
        }
        res
    }

    pub fn pow(&self, n: ModuloFieldElement) -> Self {
        assert!(n.value.is_positive());
        let mut cur_pow = self.clone();
        let mut res = ModuloFieldElement::one(Some(self.modulus.clone()));
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

// Implement conversion from BigInt for ModuloFieldElement.
impl From<BigInt> for ModuloFieldElement {
    fn from(value: BigInt) -> Self {
        ModuloFieldElement::new(value, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::ToBigInt;

    #[test]
    fn test_new_default_modulus() {
        let fe = ModuloFieldElement::new(10.to_bigint().unwrap(), None);
        assert_eq!(fe.value, 10.to_bigint().unwrap());
        assert_eq!(fe.modulus, (3_i64 * (1 << 30) + 1).to_bigint().unwrap());
    }

    #[test]
    fn test_new_custom_modulus() {
        let fe = ModuloFieldElement::new(10.to_bigint().unwrap(), Some(7.to_bigint().unwrap()));
        assert_eq!(fe.value, 10.to_bigint().unwrap() % 7.to_bigint().unwrap());
        assert_eq!(fe.modulus, 7.to_bigint().unwrap());
    }

    #[test]
    fn test_zero() {
        let fe_zero = ModuloFieldElement::zero(None);
        assert_eq!(fe_zero.value, BigInt::zero());
        assert_eq!(
            fe_zero.modulus,
            (3_i64 * (1 << 30) + 1).to_bigint().unwrap()
        );
    }

    #[test]
    fn test_one() {
        let fe_one = ModuloFieldElement::one(None);
        assert_eq!(fe_one.value, BigInt::one());
        assert_eq!(fe_one.modulus, (3_i64 * (1 << 30) + 1).to_bigint().unwrap());
    }

    #[test]
    fn test_addition() {
        let fe1 = ModuloFieldElement::new(10.to_bigint().unwrap(), None);
        let fe2 = ModuloFieldElement::new(15.to_bigint().unwrap(), None);
        let fe3 = fe1 + fe2;
        assert_eq!(fe3.value, (10 + 15).to_bigint().unwrap() % fe3.modulus);
    }

    #[test]
    fn test_subtraction() {
        let fe1 = ModuloFieldElement::new(20.to_bigint().unwrap(), None);
        let fe2 = ModuloFieldElement::new(15.to_bigint().unwrap(), None);
        let fe3 = fe1 - fe2;
        assert_eq!(fe3.value, (20 - 15).to_bigint().unwrap() % fe3.modulus);
    }

    #[test]
    fn test_multiplication() {
        let fe1 = ModuloFieldElement::new(5.to_bigint().unwrap(), None);
        let fe2 = ModuloFieldElement::new(4.to_bigint().unwrap(), None);
        let fe3 = fe1 * fe2;
        assert_eq!(fe3.value, (5 * 4).to_bigint().unwrap() % fe3.modulus);
    }

    #[test]
    fn test_inverse() {
        let fe1 = ModuloFieldElement::new(7.to_bigint().unwrap(), Some(17.to_bigint().unwrap()));
        let fe_inv = fe1.inverse();
        assert_eq!(fe_inv.value, 5.to_bigint().unwrap());
        assert_eq!((fe1 * fe_inv).value, 1.to_bigint().unwrap());
    }

    #[test]
    fn test_negation() {
        let fe1 = ModuloFieldElement::new(10.to_bigint().unwrap(), None);
        let fe_neg = -fe1;
        assert_eq!(
            fe_neg.value,
            ((-10 + DEFAULT_K_MODULES).to_bigint().unwrap()) % fe_neg.modulus
        );
    }

    #[test]
    fn test_division() {
        let fe1 = ModuloFieldElement::new(12.to_bigint().unwrap(), Some(17.to_bigint().unwrap()));
        let fe2 = ModuloFieldElement::new(3.to_bigint().unwrap(), Some(17.to_bigint().unwrap()));
        let fe_div = fe1 / fe2.clone();
        assert_eq!(
            fe_div.value,
            (12.to_bigint().unwrap() * fe2.inverse().value) % fe_div.modulus
        );
    }

    #[test]
    fn test_exponentiation() {
        let fe1 = ModuloFieldElement::new(3.to_bigint().unwrap(), Some(17.to_bigint().unwrap()));
        let fe_exp = fe1.pow_scalar(4);
        assert_eq!(
            fe_exp.value,
            3.to_bigint().unwrap().pow(4) % 17.to_bigint().unwrap()
        );
    }

    #[test]
    fn test_is_order() {
        let fe1 = ModuloFieldElement::new(2.to_bigint().unwrap(), Some(17.to_bigint().unwrap()));
        assert!(fe1.is_order(8));
        assert!(!fe1.is_order(3));
    }

    #[test]
    fn test_random_element() {
        let excluded_elements = vec![
            ModuloFieldElement::new(2.to_bigint().unwrap(), None),
            ModuloFieldElement::new(3.to_bigint().unwrap(), None),
        ];
        let fe_random = ModuloFieldElement::random_element(&excluded_elements);
        assert!(!excluded_elements.contains(&fe_random));
    }
}
