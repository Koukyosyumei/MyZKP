use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;
use std::fmt;
use std::fmt::Debug;
use std::hash::Hash;
use std::hash::Hasher;
use std::marker::PhantomData;
use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::modules::field::{Field, FiniteFieldElement, ModulusValue};
use crate::modules::polynomial::Polynomial;
use crate::modules::ring::Ring;
use crate::modules::utils::extended_euclidean;

pub trait IrreduciblePoly<F: Field>: Debug + Clone + Hash {
    fn modulus() -> Polynomial<F>;
}

#[derive(Clone, Debug)]
pub struct ExtendedFieldElement<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> {
    poly: Polynomial<FiniteFieldElement<M>>,
    _phantom: PhantomData<P>,
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> ExtendedFieldElement<M, P> {
    pub fn new(poly: Polynomial<FiniteFieldElement<M>>) -> Self {
        let mut result = Self {
            poly: poly,
            _phantom: PhantomData,
        };
        result.reduce();
        result
    }

    fn reduce(&mut self) {
        self.poly = self.poly.clone() % P::modulus();
    }

    pub fn degree(&self) -> isize {
        P::modulus().degree()
    }

    pub fn from_base_field(value: FiniteFieldElement<M>) -> Self {
        Self::new(Polynomial { coef: vec![value] })
    }
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> Hash
    for ExtendedFieldElement<M, P>
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        for v in &self.poly.coef.clone() {
            v.hash(state);
        }
    }
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> PartialEq
    for ExtendedFieldElement<M, P>
{
    fn eq(&self, other: &Self) -> bool {
        self.poly == other.poly
    }
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> Eq for ExtendedFieldElement<M, P> {}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> Zero
    for ExtendedFieldElement<M, P>
{
    fn zero() -> Self {
        ExtendedFieldElement::<M, P>::new(Polynomial::<FiniteFieldElement<M>>::zero())
    }

    fn is_zero(&self) -> bool {
        self.poly.is_zero()
    }
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> One
    for ExtendedFieldElement<M, P>
{
    fn one() -> Self {
        ExtendedFieldElement::<M, P>::new(Polynomial::<FiniteFieldElement<M>>::one())
    }
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> Add
    for ExtendedFieldElement<M, P>
{
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self::new(self.poly + other.poly)
    }
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> Sub
    for ExtendedFieldElement<M, P>
{
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self::new(self.poly - other.poly)
    }
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> Mul
    for ExtendedFieldElement<M, P>
{
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Self::new(self.poly * other.poly)
    }
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> Div
    for ExtendedFieldElement<M, P>
{
    type Output = Self;

    fn div(self, other: Self) -> Self {
        Self::new(self.poly / other.poly)
    }
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> Neg
    for ExtendedFieldElement<M, P>
{
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(-self.poly)
    }
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> Field
    for ExtendedFieldElement<M, P>
{
    fn inverse(&self) -> Self {
        // Implementation of inverse using extended Euclidean algorithm for polynomials
        let (r, s, t) = extended_euclidean(self.poly.clone(), P::modulus());
        ExtendedFieldElement::<M, P>::new(t.clone())
        //assert_eq!(gcd.degree(), 0, "Element is not invertible");
        //Self::new(s * gcd.nth_coefficient(0).inverse())
    }
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> Ring
    for ExtendedFieldElement<M, P>
{
    fn pow(&self, n: BigInt) -> Self {
        let mut result = Self::new(Polynomial::from_monomials(
            &[FiniteFieldElement::<M>::one()],
        ));
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

    fn get_value(&self) -> BigInt {
        unimplemented!("Not applicable for extended field elements")
    }

    fn from_value<V: Into<BigInt>>(_value: V) -> Self {
        unimplemented!("Use from_base_field for extended field elements")
    }

    fn random_element(_exclude_elements: &[Self]) -> Self {
        unimplemented!("Random element generation for extended field not implemented")
    }
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> fmt::Display
    for ExtendedFieldElement<M, P>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.poly)
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
        #[derive(Debug, Clone, PartialEq, Hash)]
        pub struct Ip7;

        impl IrreduciblePoly<FiniteFieldElement<Mod7>> for Ip7 {
            // x^2 + 1
            fn modulus() -> Polynomial<FiniteFieldElement<Mod7>> {
                Polynomial {
                    coef: vec![
                        FiniteFieldElement::<Mod7>::from_value(1),
                        FiniteFieldElement::<Mod7>::zero(),
                        FiniteFieldElement::<Mod7>::from_value(1),
                    ],
                }
            }
        }

        // 2x + 3
        let a = ExtendedFieldElement::<Mod7, Ip7>::new(Polynomial {
            coef: vec![
                FiniteFieldElement::from_value(3),
                FiniteFieldElement::from_value(2),
            ],
        });

        // 4x + 1
        let b = ExtendedFieldElement::<Mod7, Ip7>::new(Polynomial {
            coef: vec![
                FiniteFieldElement::from_value(1),
                FiniteFieldElement::from_value(4),
            ],
        });

        // Addition
        let sum = a.clone() + b.clone();
        assert_eq!(
            sum,
            ExtendedFieldElement::<Mod7, Ip7>::new(Polynomial {
                coef: vec![
                    FiniteFieldElement::from_value(4),
                    FiniteFieldElement::from_value(6)
                ],
            })
        );

        // Multiplication
        let product = a.clone() * b.clone();
        assert_eq!(
            product,
            ExtendedFieldElement::<Mod7, Ip7>::new(Polynomial {
                coef: vec![FiniteFieldElement::from_value(2),],
            },)
        );

        // Inverse
        let inv_a = a.clone().inverse();
        let product = a.clone() * inv_a;
        assert_eq!(
            product,
            ExtendedFieldElement::<Mod7, Ip7>::from_base_field(FiniteFieldElement::one(),)
        );

        /*
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

                },

            )
        );
                */
    }
}
