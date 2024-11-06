use crate::modules::ring::Ring;
use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;
use std::fmt;
use std::hash::Hash;
use std::hash::Hasher;
use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::modules::field::{Field, FiniteFieldElement, ModulusValue};
use crate::modules::polynomial::Polynomial;

// Assuming the FiniteFieldElement, ModulusValue, Field, and Polynomial structs are already defined

#[derive(Clone, Debug)]
pub struct ExtendedFieldElement<M: ModulusValue> {
    poly: Polynomial<FiniteFieldElement<M>>,
    irreducible_poly: Polynomial<FiniteFieldElement<M>>,
}

impl<M: ModulusValue> ExtendedFieldElement<M> {
    pub fn new(
        poly: Polynomial<FiniteFieldElement<M>>,
        irreducible_poly: Polynomial<FiniteFieldElement<M>>,
    ) -> Self {
        let mut result = Self {
            poly,
            irreducible_poly,
        };
        result.reduce();
        result
    }

    fn reduce(&mut self) {
        self.poly = self.poly.clone() % self.irreducible_poly.clone();
    }

    pub fn degree(&self) -> isize {
        self.irreducible_poly.degree()
    }

    pub fn from_base_field(
        value: FiniteFieldElement<M>,
        irreducible_poly: Polynomial<FiniteFieldElement<M>>,
    ) -> Self {
        Self::new(
            Polynomial {
                poly: vec![value],
                var: "x".to_string(),
            },
            irreducible_poly,
        )
    }
}

impl<M: ModulusValue> Hash for ExtendedFieldElement<M> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        for v in &self.poly.poly.clone() {
            v.hash(state);
        }
    }
}

impl<M: ModulusValue> PartialEq for ExtendedFieldElement<M> {
    fn eq(&self, other: &Self) -> bool {
        self.poly == other.poly
    }
}

impl<M: ModulusValue> Eq for ExtendedFieldElement<M> {}

impl<M: ModulusValue> Ring for ExtendedFieldElement<M> {
    fn zero() -> Self {
        ExtendedFieldElement::<M>::new(Polynomial::<FiniteFieldElement<M>>::zero())
    }

    fn one() -> Self {
        ExtendedFieldElement::<M>::new(Polynomial::<FiniteFieldElement<M>>::one())
    }
}

impl<M: ModulusValue> Mul<i64> for ExtendedFieldElement<M> {
    type Output = Self;

    fn mul(self, n: i64) -> Self {
        let modulus = M::modulus();
        ExtendedFieldElement::<M>::new((&self.value * n.to_bigint().unwrap()) % &modulus)
    }
}

impl<M: ModulusValue> Add for ExtendedFieldElement<M> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        assert_eq!(
            self.irreducible_poly, other.irreducible_poly,
            "Incompatible field extensions"
        );
        Self::new(self.poly + other.poly, self.irreducible_poly)
    }
}

impl<M: ModulusValue> Sub for ExtendedFieldElement<M> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        assert_eq!(
            self.irreducible_poly, other.irreducible_poly,
            "Incompatible field extensions"
        );
        Self::new(self.poly - other.poly, self.irreducible_poly)
    }
}

impl<M: ModulusValue> Mul for ExtendedFieldElement<M> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        assert_eq!(
            self.irreducible_poly, other.irreducible_poly,
            "Incompatible field extensions"
        );
        Self::new(self.poly * other.poly, self.irreducible_poly)
    }
}

impl<M: ModulusValue> Div for ExtendedFieldElement<M> {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        assert_eq!(
            self.irreducible_poly, other.irreducible_poly,
            "Incompatible field extensions"
        );
        Self::new(self.poly / other.poly, self.irreducible_poly)
    }
}

impl<M: ModulusValue> Neg for ExtendedFieldElement<M> {
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(-self.poly, self.irreducible_poly)
    }
}

impl<M: ModulusValue> Field for ExtendedFieldElement<M> {
    fn inverse(&self) -> Self {
        // Implementation of inverse using extended Euclidean algorithm for polynomials
        let (gcd, s, _) = extended_euclidean(&self.poly, &self.irreducible_poly);
        assert_eq!(gcd.degree(), 0, "Element is not invertible");
        Self::new(
            s * gcd.nth_coefficient(0).inverse(),
            self.irreducible_poly.clone(),
        )
    }

    fn get_value(&self) -> BigInt {
        unimplemented!("Not applicable for extended field elements")
    }

    fn pow(&self, n: BigInt) -> Self {
        let mut result = Self::new(
            Polynomial::from_monomials(&[FiniteFieldElement::<M>::one()]),
            self.irreducible_poly.clone(),
        );
        let mut base = self.clone();
        let mut exp = n;

        while exp > BigInt::zero() {
            if &exp % BigInt::from(2) == BigInt::one() {
                result = result * base.clone();
            }
            base = base.clone() * base;
            exp /= 2;
        }

        result
    }

    fn from_value<V: Into<BigInt>>(value: V) -> Self {
        unimplemented!("Use from_base_field for extended field elements")
    }

    fn random_element(exclude_elements: &[Self]) -> Self {
        unimplemented!("Random element generation for extended field not implemented")
    }
}

impl<M: ModulusValue> fmt::Display for ExtendedFieldElement<M> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.poly)
    }
}

// Helper function for extended Euclidean algorithm for polynomials
fn extended_euclidean<F: Field>(
    a: &Polynomial<F>,
    b: &Polynomial<F>,
) -> (Polynomial<F>, Polynomial<F>, Polynomial<F>) {
    if b.degree() == -1 {
        (
            a.clone(),
            Polynomial::from_monomials(&[F::one()]),
            Polynomial::from_monomials(&[F::zero()]),
        )
    } else {
        let (q, r) = (a.clone() / b.clone(), a.clone() % b.clone());
        let (gcd, s, t) = extended_euclidean(b, &r);
        (gcd, t.clone(), s - q * t)
    }
}

// Test the implementation
#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    use crate::define_myzkp_modulus_type;

    define_myzkp_modulus_type!(Mod7, "7");

    #[test]
    fn test_extended_field_operations() {
        // Define the irreducible polynomial x^2 + 1
        let irreducible_poly = Polynomial {
            poly: vec![
                FiniteFieldElement::<Mod7>::from_value(1),
                FiniteFieldElement::<Mod7>::zero(),
                FiniteFieldElement::<Mod7>::from_value(1),
            ],
            var: "x".to_string(),
        };

        let a = ExtendedFieldElement::<Mod7>::new(
            Polynomial {
                poly: vec![
                    FiniteFieldElement::from_value(3),
                    FiniteFieldElement::from_value(2),
                ],
                var: "x".to_string(),
            },
            irreducible_poly.clone(),
        );
        let b = ExtendedFieldElement::<Mod7>::new(
            Polynomial {
                poly: vec![
                    FiniteFieldElement::from_value(1),
                    FiniteFieldElement::from_value(4),
                ],
                var: "x".to_string(),
            },
            irreducible_poly.clone(),
        );

        // Addition
        let sum = a.clone() + b.clone();
        assert_eq!(
            sum,
            ExtendedFieldElement::<Mod7>::new(
                Polynomial {
                    poly: vec![
                        FiniteFieldElement::from_value(4),
                        FiniteFieldElement::from_value(6)
                    ],
                    var: "x".to_string(),
                },
                irreducible_poly.clone()
            )
        );

        // Multiplication
        let product = a.clone() * b.clone();
        assert_eq!(
            product,
            ExtendedFieldElement::<Mod7>::new(
                Polynomial {
                    poly: vec![
                        FiniteFieldElement::from_value(2),
                        FiniteFieldElement::from_value(5)
                    ],
                    var: "x".to_string(),
                },
                irreducible_poly.clone()
            )
        );

        // Inverse
        let inv_a = a.inverse();
        let product = a * inv_a;
        assert_eq!(
            product,
            ExtendedFieldElement::<Mod7>::from_base_field(
                FiniteFieldElement::one(),
                irreducible_poly.clone()
            )
        );

        // Exponentiation
        let a_cubed = a.pow(BigInt::from(3));
        assert_eq!(
            a_cubed,
            ExtendedFieldElement::<Mod7>::new(
                Polynomial {
                    poly: vec![
                        FiniteFieldElement::from_value(5),
                        FiniteFieldElement::from_value(6)
                    ],
                    var: "x".to_string(),
                },
                irreducible_poly.clone()
            )
        );
    }
}
