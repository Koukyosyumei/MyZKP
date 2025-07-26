//! # Extended Finite Field Element Implementation
//!
//! This module provides an implementation of extended finite field elements,
//! which are built on top of base finite fields. It includes structures and
//! traits for defining and working with extended finite fields, as well as
//! implementations of arithmetic operations in these fields.
//!
//! ## Key Components
//!
//! - `IrreduciblePoly` trait: Defines the irreducible polynomial for the field extension.
//! - `ExtendedFieldElement<M, P>` struct: Represents an element in an extended finite field.
//! - Arithmetic operations: Addition, subtraction, multiplication, and division.
//!
//! ## Features
//!
//! - Creation of extended field elements from base field polynomials.
//! - Arithmetic operations in the extended field.
//! - Reduction of elements modulo the irreducible polynomial.
//! - Inverse computation in the extended field.
//! - Conversion between base field elements and extended field elements.
//!
//! ## Usage
//!
//! To use this module, define a base field, an irreducible polynomial, and then
//! create `ExtendedFieldElement` instances:
//!
//! ```
//! use std::str::FromStr;
//! use paste::paste;
//! use num_bigint::BigInt;
//! use lazy_static::lazy_static;
//! use serde::{Deserialize, Serialize};
//! use myzkp::define_myzkp_modulus_type;
//! use myzkp::define_extension_field;
//! use myzkp::modules::algebra::ring::Ring;
//! use myzkp::modules::algebra::field::ModulusValue;
//! use myzkp::modules::algebra::field::FiniteFieldElement;
//! use myzkp::modules::algebra::polynomial::Polynomial;
//! use myzkp::modules::algebra::efield::IrreduciblePoly;
//! use myzkp::modules::algebra::efield::ExtendedFieldElement;
//!
//! define_myzkp_modulus_type!(Mod7, "7");
//! define_extension_field!(
//!     Ip7,
//!     FiniteFieldElement<Mod7>,
//!     Polynomial {
//!         coef: vec![
//!             FiniteFieldElement::<Mod7>::from_value(1),
//!             FiniteFieldElement::<Mod7>::from_value(0),
//!             FiniteFieldElement::<Mod7>::from_value(1),
//!         ],
//!     }
//! );
//!
//! let a = ExtendedFieldElement::<Mod7, Ip7>::new(Polynomial {
//!     coef: vec![
//!         FiniteFieldElement::from_value(2),
//!         FiniteFieldElement::from_value(1),
//!     ],
//! });
//! ```
//!
//! ## Note
//!
//! This implementation builds upon the base finite field implementation and uses
//! the `Polynomial` struct for representing elements. It is designed to work with
//! various base fields and irreducible polynomials, allowing for flexible creation
//! of field extensions.

use std::fmt;
use std::fmt::Debug;
use std::hash::Hash;
use std::hash::Hasher;
use std::marker::PhantomData;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;
use serde::{Deserialize, Serialize};

use crate::modules::algebra::field::{Field, FiniteFieldElement, ModulusValue};
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::ring::Ring;

pub trait IrreduciblePoly<F: Field>: Debug + Clone + Hash {
    fn modulus() -> &'static Polynomial<F>;
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ExtendedFieldElement<M: ModulusValue, P: IrreduciblePoly<FiniteFieldElement<M>>> {
    pub poly: Polynomial<FiniteFieldElement<M>>,
    _phantom: PhantomData<P>,
}

impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>>
    ExtendedFieldElement<M, P>
{
    pub fn new(poly: Polynomial<FiniteFieldElement<M>>) -> Self {
        let result = Self {
            poly: poly,
            _phantom: PhantomData,
        };
        result.sanitize()
    }

    pub fn degree(&self) -> isize {
        P::modulus().degree()
    }

    pub fn from_base_field(value: FiniteFieldElement<M>) -> Self {
        Self::new((Polynomial { coef: vec![value] }).reduce()).sanitize()
    }
}

impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>> Field
    for ExtendedFieldElement<M, P>
{
    fn inverse(&self) -> Self {
        //if self.poly.is_zero() {
        //    return None; // Zero has no inverse
        //}

        let mut lm = Polynomial::<FiniteFieldElement<M>>::one();
        let mut hm = Polynomial::zero();
        let mut low = self.poly.clone();
        let mut high = P::modulus().clone();

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

    fn add_m1_ref(&self, _other: &Self) -> Self {
        unimplemented!("Not applicable for extended field elements")
    }

    fn mul_m1_ref(&self, _other: &Self) -> Self {
        unimplemented!("Not applicable for extended field elements")
    }

    fn sub_m1_ref(&self, _other: &Self) -> Self {
        unimplemented!("Not applicable for extended field elements")
    }

    fn pow_m1<V: Into<BigInt>>(&self, _n: V) -> Self {
        unimplemented!("Not applicable for extended field elements")
    }

    fn sanitize(&self) -> Self {
        Self {
            poly: &self.poly.reduce() % P::modulus(),
            _phantom: PhantomData,
        }
    }

    fn sample(byte_array: &[u8]) -> Self {
        let mut acc: usize = 0;
        for &b in byte_array {
            acc = (acc << 8) ^ (b as usize);
        }
        Self::from_value(acc)
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

impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>> PartialEq
    for ExtendedFieldElement<M, P>
{
    fn eq(&self, other: &Self) -> bool {
        self.poly == other.poly
    }
}

impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>> Eq
    for ExtendedFieldElement<M, P>
{
}

impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>> Zero
    for ExtendedFieldElement<M, P>
{
    fn zero() -> Self {
        ExtendedFieldElement::<M, P>::new(Polynomial::<FiniteFieldElement<M>>::zero())
    }

    fn is_zero(&self) -> bool {
        self.poly.is_zero()
    }
}

impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>> One
    for ExtendedFieldElement<M, P>
{
    fn one() -> Self {
        ExtendedFieldElement::<M, P>::new(Polynomial::<FiniteFieldElement<M>>::one())
    }
}

impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>> Add
    for ExtendedFieldElement<M, P>
{
    type Output = Self;

    fn add(self, other: Self) -> Self {
        self.add_ref(&other)
    }
}

impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>> AddAssign
    for ExtendedFieldElement<M, P>
{
    fn add_assign(&mut self, other: Self) {
        self.add_assign_ref(&other)
    }
}

impl<'a, M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>>
    Add<&'a ExtendedFieldElement<M, P>> for ExtendedFieldElement<M, P>
{
    type Output = ExtendedFieldElement<M, P>;

    fn add(self, other: &'a ExtendedFieldElement<M, P>) -> ExtendedFieldElement<M, P> {
        self.add_ref(other)
    }
}

impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>> Sub
    for ExtendedFieldElement<M, P>
{
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self.sub_ref(&other)
    }
}

impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>> SubAssign
    for ExtendedFieldElement<M, P>
{
    fn sub_assign(&mut self, other: Self) {
        self.sub_assign_ref(&other)
    }
}

impl<'a, M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>>
    Sub<&'a ExtendedFieldElement<M, P>> for ExtendedFieldElement<M, P>
{
    type Output = ExtendedFieldElement<M, P>;

    fn sub(self, other: &'a ExtendedFieldElement<M, P>) -> ExtendedFieldElement<M, P> {
        self.sub_ref(other)
    }
}

impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>> Mul
    for ExtendedFieldElement<M, P>
{
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        self.mul_ref(&other)
    }
}

impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>> MulAssign
    for ExtendedFieldElement<M, P>
{
    fn mul_assign(&mut self, other: Self) {
        self.mul_assign_ref(&other)
    }
}

impl<'a, M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>>
    Mul<&'a ExtendedFieldElement<M, P>> for ExtendedFieldElement<M, P>
{
    type Output = ExtendedFieldElement<M, P>;

    fn mul(self, other: &'a ExtendedFieldElement<M, P>) -> ExtendedFieldElement<M, P> {
        self.mul_ref(other)
    }
}

impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>> Div
    for ExtendedFieldElement<M, P>
{
    type Output = Self;

    fn div(self, other: Self) -> Self {
        self.div_ref(&other)
    }
}

impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>> Neg
    for ExtendedFieldElement<M, P>
{
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(-self.poly)
    }
}

impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>> Ring
    for ExtendedFieldElement<M, P>
{
    fn add_ref(&self, other: &Self) -> Self {
        Self::new(&self.poly + &other.poly)
    }

    fn add_assign_ref(&mut self, other: &Self) {
        self.poly += &other.poly;
        self.sanitize();
    }

    fn mul_ref(&self, other: &Self) -> Self {
        Self::new(&self.poly * &other.poly)
    }

    fn mul_assign_ref(&mut self, other: &Self) {
        self.poly *= &other.poly;
        self.sanitize();
    }

    fn sub_ref(&self, other: &Self) -> Self {
        Self::new(&self.poly - &other.poly)
    }

    fn sub_assign_ref(&mut self, other: &Self) {
        self.poly -= &other.poly;
        self.sanitize();
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

impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>> fmt::Display
    for ExtendedFieldElement<M, P>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.poly)
    }
}

#[macro_export]
macro_rules! define_extension_field {
    ($name:ident, $base_field:ty, $modulus:expr) => {
        paste! {#[derive(Debug, Clone, PartialEq, Hash)]
            pub struct $name;

            lazy_static! {
                static ref [<MODULUS_ $name>]: Polynomial<$base_field> = $modulus;
            }

            impl IrreduciblePoly<$base_field> for $name {
                fn modulus() -> &'static Polynomial<$base_field> {
                    &[<MODULUS_ $name>]
                }
            }
        }
    };
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

    define_extension_field!(
        Ip2,
        FiniteFieldElement<Mod2>,
        Polynomial {
            coef: vec![
                FiniteFieldElement::<Mod2>::from_value(1),
                FiniteFieldElement::<Mod2>::from_value(1),
                FiniteFieldElement::<Mod2>::from_value(1),
            ],
        }
    );

    define_extension_field!(
        Ip7,
        FiniteFieldElement<Mod7>,
        Polynomial {
            coef: vec![
                FiniteFieldElement::<Mod7>::from_value(1),
                FiniteFieldElement::<Mod7>::from_value(0),
                FiniteFieldElement::<Mod7>::from_value(1),
            ],
        }
    );

    #[test]
    fn test_extended_field_operations_mod2() {
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
