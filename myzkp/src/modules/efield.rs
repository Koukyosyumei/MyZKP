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
    pub poly: Polynomial<FiniteFieldElement<M>>,
    _phantom: PhantomData<P>,
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> ExtendedFieldElement<M, P> {
    pub fn new(poly: Polynomial<FiniteFieldElement<M>>) -> Self {
        let result = Self {
            poly: poly,
            _phantom: PhantomData,
        };
        result.reduce()
    }

    fn reduce(&self) -> Self {
        Self {
            poly: self.poly.reduce() % P::modulus(),
            _phantom: PhantomData,
        }
    }

    pub fn degree(&self) -> isize {
        P::modulus().degree()
    }

    pub fn from_base_field(value: FiniteFieldElement<M>) -> Self {
        Self::new((Polynomial { coef: vec![value] }).reduce()).reduce()
    }
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> Field
    for ExtendedFieldElement<M, P>
{
    fn inverse(&self) -> Self {
        //if self.poly.is_zero() {
        //    return None; // Zero has no inverse
        //}

        let mut lm = Polynomial::<FiniteFieldElement<M>>::one();
        let mut hm = Polynomial::zero();
        let mut low = self.poly.clone();
        let mut high = P::modulus();

        while !low.is_zero() {
            let q = &high / &low;
            let r = &high % &low;
            let nm = hm - (&lm * &q);
            high = low;
            hm = lm;
            low = r;
            lm = nm;
        }

        //if high.degree() != 0 {
        //    return None; // Not invertible
        //}

        Self::new(hm * high.coef[0].inverse())
    }

    fn div_ref(&self, other: &Self) -> Self {
        self.mul_ref(&other.inverse())
    }
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> Hash
    for ExtendedFieldElement<M, P>
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        for v in &self.poly.coef {
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
        self.add_ref(&other)
    }
}

impl<'a, M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>>
    Add<&'a ExtendedFieldElement<M, P>> for ExtendedFieldElement<M, P>
{
    type Output = ExtendedFieldElement<M, P>;

    fn add(self, other: &'a ExtendedFieldElement<M, P>) -> ExtendedFieldElement<M, P> {
        self.add_ref(other)
    }
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> Sub
    for ExtendedFieldElement<M, P>
{
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self.sub_ref(&other)
    }
}

impl<'a, M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>>
    Sub<&'a ExtendedFieldElement<M, P>> for ExtendedFieldElement<M, P>
{
    type Output = ExtendedFieldElement<M, P>;

    fn sub(self, other: &'a ExtendedFieldElement<M, P>) -> ExtendedFieldElement<M, P> {
        self.sub_ref(other)
    }
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> Mul
    for ExtendedFieldElement<M, P>
{
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        self.mul_ref(&other)
    }
}

impl<'a, M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>>
    Mul<&'a ExtendedFieldElement<M, P>> for ExtendedFieldElement<M, P>
{
    type Output = ExtendedFieldElement<M, P>;

    fn mul(self, other: &'a ExtendedFieldElement<M, P>) -> ExtendedFieldElement<M, P> {
        self.mul_ref(other)
    }
}

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> Div
    for ExtendedFieldElement<M, P>
{
    type Output = Self;

    fn div(self, other: Self) -> Self {
        self.div_ref(&other)
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

impl<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> Ring
    for ExtendedFieldElement<M, P>
{
    fn add_ref(&self, other: &Self) -> Self {
        Self::new(&self.poly + &other.poly)
    }

    fn mul_ref(&self, other: &Self) -> Self {
        Self::new(&self.poly * &other.poly)
    }

    fn sub_ref(&self, other: &Self) -> Self {
        Self::new(&self.poly - &other.poly)
    }

    fn pow<V: Into<BigInt>>(&self, n: V) -> Self {
        let mut base = self.clone();
        let mut exponent: BigInt = n.into();

        let mut result = Self::one();
        while exponent > BigInt::zero() {
            if &exponent % BigInt::from(2) == BigInt::one() {
                result = result.mul_ref(&base);
            }
            exponent /= 2;
            base = base.mul_ref(&base);
        }

        result
    }

    fn get_value(&self) -> BigInt {
        unimplemented!("Not applicable for extended field elements")
    }

    fn from_value<V: Into<BigInt>>(value: V) -> Self {
        ExtendedFieldElement::<M, P>::new(
            Polynomial::<FiniteFieldElement<M>>::one() * FiniteFieldElement::<M>::from_value(value),
        )
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
    use lazy_static::lazy_static;
    use paste::paste;
    use std::str::FromStr;

    use crate::define_myzkp_modulus_type;

    define_myzkp_modulus_type!(Mod2, "2");
    define_myzkp_modulus_type!(Mod7, "7");

    #[test]
    fn test_extended_field_operations_mod2() {
        #[derive(Debug, Clone, PartialEq, Hash)]
        pub struct Ip2;

        impl IrreduciblePoly<FiniteFieldElement<Mod2>> for Ip2 {
            // x^2 + x + 1
            fn modulus() -> Polynomial<FiniteFieldElement<Mod2>> {
                Polynomial {
                    coef: vec![
                        FiniteFieldElement::<Mod2>::from_value(1),
                        FiniteFieldElement::<Mod2>::from_value(1),
                        FiniteFieldElement::<Mod2>::from_value(1),
                    ],
                }
            }
        }

        // x + 1
        let a = ExtendedFieldElement::<Mod2, Ip2>::new(Polynomial {
            coef: vec![
                FiniteFieldElement::from_value(1),
                FiniteFieldElement::from_value(1),
            ],
        });

        // x
        let b = ExtendedFieldElement::<Mod2, Ip2>::new(Polynomial {
            coef: vec![
                FiniteFieldElement::from_value(0),
                FiniteFieldElement::from_value(1),
            ],
        });

        // Addition
        let sum = a.add_ref(&b);
        assert_eq!(
            sum,
            ExtendedFieldElement::<Mod2, Ip2>::new(Polynomial {
                coef: vec![FiniteFieldElement::from_value(1)],
            })
        );

        // Multiplication
        let product = a.mul_ref(&b);
        assert_eq!(
            product,
            ExtendedFieldElement::<Mod2, Ip2>::new(Polynomial {
                coef: vec![FiniteFieldElement::from_value(1)],
            },)
        );

        // Inverse
        let inv_a = a.inverse();
        let product = a.mul_ref(&inv_a);
        assert_eq!(
            product,
            ExtendedFieldElement::<Mod2, Ip2>::from_base_field(FiniteFieldElement::one())
        );

        assert_eq!(
            a.div_ref(&a),
            ExtendedFieldElement::<Mod2, Ip2>::from_base_field(FiniteFieldElement::one())
        );
    }

    #[test]
    fn test_extended_field_operations_mod7() {
        #[derive(Debug, Clone, PartialEq, Hash)]
        pub struct Ip7;

        impl IrreduciblePoly<FiniteFieldElement<Mod7>> for Ip7 {
            // x^2 + x + 1
            fn modulus() -> Polynomial<FiniteFieldElement<Mod7>> {
                Polynomial {
                    coef: vec![
                        FiniteFieldElement::<Mod7>::from_value(1),
                        FiniteFieldElement::<Mod7>::from_value(0),
                        FiniteFieldElement::<Mod7>::from_value(1),
                    ],
                }
            }
        }

        // x + 2
        let a = ExtendedFieldElement::<Mod7, Ip7>::new(Polynomial {
            coef: vec![
                FiniteFieldElement::from_value(2),
                FiniteFieldElement::from_value(1),
            ],
        });

        // 4x + 6
        let b = ExtendedFieldElement::<Mod7, Ip7>::new(Polynomial {
            coef: vec![
                FiniteFieldElement::from_value(6),
                FiniteFieldElement::from_value(4),
            ],
        });

        // Inverse
        let inv_a = a.inverse();
        assert_eq!(inv_a, b);

        let product = a.mul_ref(&inv_a);
        assert_eq!(
            product,
            ExtendedFieldElement::<Mod7, Ip7>::from_base_field(FiniteFieldElement::one())
        );
    }
}
