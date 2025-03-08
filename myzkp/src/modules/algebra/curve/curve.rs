use std::fmt;
use std::fmt::Debug;
use std::marker::PhantomData;
use std::ops::{Add, AddAssign, Mul, Neg, Sub};

use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

use crate::modules::algebra::field::{Field, FiniteFieldElement, ModulusValue};

pub trait EllipticCurve: Debug + Clone + PartialEq {
    fn get_a() -> BigInt;
    fn get_b() -> BigInt;
}

#[derive(Debug, Clone, PartialEq)]
pub struct EllipticCurvePoint<F: Field, E: EllipticCurve> {
    pub x: Option<F>,
    pub y: Option<F>,
    _phantom: PhantomData<E>,
}

impl<F: Field, E: EllipticCurve> EllipticCurvePoint<F, E> {
    pub fn new(x: F, y: F) -> Self {
        // let a = E::get_a();
        // let b = E::get_b();
        // assert!(y.pow(2) == x.pow(3) + a * x + b, "Point is not on the curve");
        EllipticCurvePoint {
            x: Some(x),
            y: Some(y),
            _phantom: PhantomData,
        }
    }

    pub fn point_at_infinity() -> Self {
        EllipticCurvePoint {
            x: None,
            y: None,
            _phantom: PhantomData,
        }
    }

    pub fn is_point_at_infinity(&self) -> bool {
        self.x.is_none() || self.y.is_none()
    }

    pub fn inverse(&self) -> Self {
        if self.is_point_at_infinity() {
            return self.clone();
        }

        EllipticCurvePoint::new(self.x.clone().unwrap(), -(self.y.clone().unwrap()))
    }

    pub fn line_slope(&self, other: &Self) -> F {
        let a = F::from_value(E::get_a());

        let x1 = self.x.as_ref().unwrap();
        let y1 = self.y.as_ref().unwrap();
        let x2 = other.x.as_ref().unwrap();
        let y2 = other.y.as_ref().unwrap();

        if self.x == other.x {
            ((x1.mul_ref(&x1)) * &(F::from_value(3_i64)) + &a)
                / (y1.mul_ref(&(F::from_value(2_i64))))
        } else {
            (y2.sub_ref(&y1)) / (x2.sub_ref(&x1))
        }
    }

    pub fn double(&self) -> Self {
        if self.is_point_at_infinity() {
            return self.clone();
        }

        let slope = self.line_slope(&self);
        let x = self.x.as_ref().unwrap();
        let y = self.y.as_ref().unwrap();

        let new_x = slope.mul_ref(&slope).sub_ref(&x).sub_ref(&x);
        let new_y = -slope.mul_ref(&new_x) + slope * x - y;

        Self::new(new_x, new_y)
    }

    pub fn inplace_double(&mut self) {
        if self.is_point_at_infinity() {
            return;
        }

        let slope = self.line_slope(&self);
        let x = self.x.as_ref().unwrap();
        let y = self.y.as_ref().unwrap();

        let new_x = slope.mul_ref(&slope).sub_ref(&x).sub_ref(&x);
        let new_y = -slope.mul_ref(&new_x) + slope * x - y;

        self.x = Some(new_x);
        self.y = Some(new_y);
    }

    pub fn add_ref(&self, other: &Self) -> Self {
        if self.is_point_at_infinity() {
            return other.clone();
        }
        if other.is_point_at_infinity() {
            return self.clone();
        }

        if self.x == other.x && self.y == other.y {
            return self.double();
        } else if self.x == other.x {
            return Self::point_at_infinity();
        }

        let slope = self.line_slope(&other);
        let x1 = self.x.as_ref().unwrap();
        let y1 = self.y.as_ref().unwrap();
        let x2 = other.x.as_ref().unwrap();
        //let y2 = other.y.as_ref().unwrap();

        let new_x = slope.mul_ref(&slope).sub_ref(&x1).sub_ref(&x2);
        let new_y = ((-slope.clone()).mul_ref(&new_x)) + (&slope.mul_ref(&x1).sub_ref(&y1));
        //assert!(new_y == -slope.clone() * &new_x + slope.mul_ref(&x2).sub_ref(&y2));

        Self::new(new_x, new_y)
    }

    pub fn add_assign_ref(&mut self, other: &Self) {
        if self.is_point_at_infinity() {
            *self = other.clone();
            return;
        }
        if other.is_point_at_infinity() {
            return;
        }

        if self.x == other.x && self.y == other.y {
            *self = self.double();
            return;
        } else if self.x == other.x {
            *self = Self::point_at_infinity();
            return;
        }

        let slope = self.line_slope(other);
        let x1 = self.x.as_mut().unwrap();
        let y1 = self.y.as_mut().unwrap();
        let x2 = other.x.as_ref().unwrap();
        //let y2 = other.y.as_ref().unwrap();

        let new_x = slope.mul_ref(&slope).sub_ref(x1).sub_ref(x2);
        let new_y = (-slope.clone())
            .mul_ref(&new_x)
            .add_ref(&slope.mul_ref(x1).sub_ref(y1));
        //assert!(new_y == -slope.clone() * &new_x + slope.mul_ref(&x2).sub_ref(&y2));

        self.x = Some(new_x);
        self.y = Some(new_y);
    }

    pub fn mul_ref<V: Into<BigInt>>(&self, scalar_val: V) -> Self {
        let scalar: BigInt = scalar_val.into();
        self.mul_ref_bigint(&scalar)
    }

    pub fn mul_ref_bigint(&self, scalar: &BigInt) -> Self {
        if scalar.is_zero() {
            // Return the point at infinity for scalar * 0
            return EllipticCurvePoint::point_at_infinity();
        }

        if scalar.is_negative() {
            panic!("multiplier should be non-negative");
        }

        let mut result = EllipticCurvePoint::point_at_infinity();
        let mut current = self.clone(); // Start with the current point
        let mut scalar_bits = scalar.clone();

        while !scalar_bits.is_zero() {
            if scalar_bits.bit(0) {
                result.add_assign_ref(&current);
            }
            current.inplace_double();
            scalar_bits >>= 1; // Move to the next bit
        }

        result
    }
}

impl<F: Field, E: EllipticCurve> Add for EllipticCurvePoint<F, E> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        self.add_ref(&other)
    }
}

impl<F: Field, E: EllipticCurve> AddAssign for EllipticCurvePoint<F, E> {
    fn add_assign(&mut self, other: Self) {
        self.add_assign_ref(&other)
    }
}

impl<F: Field, E: EllipticCurve> Add for &EllipticCurvePoint<F, E> {
    type Output = EllipticCurvePoint<F, E>;

    fn add(self, other: Self) -> EllipticCurvePoint<F, E> {
        self.add_ref(&other)
    }
}

impl<F: Field, E: EllipticCurve> Neg for EllipticCurvePoint<F, E> {
    type Output = Self;
    fn neg(self) -> Self {
        if self.is_point_at_infinity() {
            return self;
        } else {
            EllipticCurvePoint::new(self.x.clone().unwrap(), -self.y.clone().unwrap())
        }
    }
}

impl<F: Field, E: EllipticCurve> Sub for EllipticCurvePoint<F, E> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self + (-other)
    }
}

impl<F: Field, E: EllipticCurve, V: Into<BigInt>> Mul<V> for EllipticCurvePoint<F, E> {
    type Output = Self;

    fn mul(self, scalar_val: V) -> Self {
        self.mul_ref(scalar_val)
    }
}

impl<F: Field, E: EllipticCurve, V: Into<BigInt>> Mul<V> for &EllipticCurvePoint<F, E> {
    type Output = EllipticCurvePoint<F, E>;

    fn mul(self, scalar_val: V) -> EllipticCurvePoint<F, E> {
        self.mul_ref(scalar_val)
    }
}

impl<F: Field, E: EllipticCurve, M: ModulusValue> Mul<FiniteFieldElement<M>>
    for &EllipticCurvePoint<F, E>
{
    type Output = EllipticCurvePoint<F, E>;

    fn mul(self, field_val: FiniteFieldElement<M>) -> EllipticCurvePoint<F, E> {
        self.mul_ref_bigint(&field_val.value)
    }
}

impl<'a, F: Field, E: EllipticCurve, M: ModulusValue> Mul<&'a FiniteFieldElement<M>>
    for &EllipticCurvePoint<F, E>
{
    type Output = EllipticCurvePoint<F, E>;

    fn mul(self, field_val: &'a FiniteFieldElement<M>) -> EllipticCurvePoint<F, E> {
        self.mul_ref_bigint(&field_val.value)
    }
}

impl<F: Field, E: EllipticCurve> fmt::Display for EllipticCurvePoint<F, E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_point_at_infinity() {
            write!(f, "(x=∞, y=∞)")
        } else {
            match (self.x.as_ref(), self.y.as_ref()) {
                (Some(x), Some(y)) => write!(f, "(x={}, y={})", x, y),
                (Some(x), None) => write!(f, "({}, ∅)", x),
                (None, Some(y)) => write!(f, "(∅, {})", y),
                (None, None) => unreachable!(), // This case is handled by is_infinity()
            }
        }
    }
}

pub fn get_lambda<F: Field, E: EllipticCurve>(
    p: &EllipticCurvePoint<F, E>,
    q: &EllipticCurvePoint<F, E>,
    r: &EllipticCurvePoint<F, E>,
) -> F {
    if p.is_point_at_infinity() || q.is_point_at_infinity() || r.is_point_at_infinity() {
        return F::one();
    }

    let p_x = p.x.as_ref().unwrap();
    let p_y = p.y.as_ref().unwrap();
    let q_x = q.x.as_ref().unwrap();
    // let q_y = q.y.clone().unwrap();
    let r_x = r.x.as_ref().unwrap();
    let r_y = r.y.as_ref().unwrap();

    if (p == q && *p_y == F::zero()) || (p != q && *p_x == *q_x) {
        return r_x.sub_ref(&p_x);
    }
    let slope = p.line_slope(&q);
    let numerator = (r_y.sub_ref(&p_y)).sub_ref(&slope.mul_ref(&(r_x.sub_ref(&p_x))));
    let denominator = r_x
        .add_ref(&p_x)
        .add_ref(&q_x)
        .sub_ref(&slope.mul_ref(&slope));
    return numerator / denominator;
}

pub fn miller<F: Field, E: EllipticCurve>(
    p: &EllipticCurvePoint<F, E>,
    q: &EllipticCurvePoint<F, E>,
    m: &BigInt,
) -> (F, EllipticCurvePoint<F, E>) {
    if p.is_point_at_infinity() || q.is_point_at_infinity() {
        return (F::one(), EllipticCurvePoint::point_at_infinity());
    }

    if p == q {
        return (F::one(), p.clone());
    }

    let mut f = F::one();
    let mut t = p.clone();

    for i in (0..(m.bits() - 1)).rev() {
        f = (f.mul_ref(&f)) * (get_lambda(&t, &t, &q));
        t = t.add_ref(&t);
        if m.bit(i) {
            f = f * (get_lambda(&t, &p, &q));
            t = t.add_ref(&p);
        }
    }

    (f, t)
}

pub fn weil_pairing<F: Field, E: EllipticCurve>(
    p: &EllipticCurvePoint<F, E>,
    q: &EllipticCurvePoint<F, E>,
    m: &BigInt,
    s: Option<&EllipticCurvePoint<F, E>>,
) -> F {
    if p.is_point_at_infinity() || q.is_point_at_infinity() {
        return F::one();
    }

    let s_value = s.unwrap();
    let (fp_qs, _) = miller(&p, &q.add_ref(&s_value), m);
    let (fp_s, _) = miller(&p, &s_value, m);
    let (fq_ps, _) = miller(&q, &p.add_ref(&(s_value.clone().neg())), m);
    let (fq_s, _) = miller(&q, &(-s_value.clone()), m);

    return (fp_qs / fp_s) / (fq_ps / fq_s);
}

pub fn general_tate_pairing<F: Field, E: EllipticCurve>(
    p: &EllipticCurvePoint<F, E>,
    q: &EllipticCurvePoint<F, E>,
    ell: &BigInt,
    modulus: &BigInt,
    s: Option<&EllipticCurvePoint<F, E>>,
) -> F {
    if p.is_point_at_infinity() || q.is_point_at_infinity() {
        return F::one();
    }

    let s_value = s.unwrap();
    let (fp_qs, _) = miller(&p, &q.add_ref(&s_value), ell);
    let (fp_s, _) = miller(&p, &s_value, ell);
    let f = fp_qs.div_ref(&fp_s);

    return f.pow((modulus - BigInt::one()) / ell);
}

pub fn tate_pairing<F: Field, E: EllipticCurve>(
    p: &EllipticCurvePoint<F, E>,
    q: &EllipticCurvePoint<F, E>,
    ell: &BigInt,
    modulus: &BigInt,
) -> F {
    if p.is_point_at_infinity() || q.is_point_at_infinity() {
        return F::one();
    }

    let (fp_q, _) = miller(p, q, ell);

    return fp_q.pow((modulus - BigInt::one()) / ell);
}

#[macro_export]
macro_rules! define_myzkp_curve_type {
    ($name:ident, $a:expr, $b:expr) => {
        #[derive(Debug, Clone, PartialEq)]
        pub struct $name;

        impl EllipticCurve for $name {
            fn get_a() -> BigInt {
                BigInt::from_str($a).unwrap()
            }
            fn get_b() -> BigInt {
                BigInt::from_str($b).unwrap()
            }
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::str::FromStr;

    use lazy_static::lazy_static;
    use num_bigint::{BigInt, ToBigInt};
    use paste::paste;
    use serde::Serialize;

    use crate::define_myzkp_modulus_type;
    use crate::modules::algebra::field::{FiniteFieldElement, ModulusValue};
    use crate::modules::algebra::ring::Ring;

    define_myzkp_modulus_type!(Mod631, "631");
    define_myzkp_curve_type!(CurveA30B34, "30", "34");

    #[test]
    fn test_weil_pairing() {
        let p = EllipticCurvePoint::<FiniteFieldElement<Mod631>, CurveA30B34>::new(
            FiniteFieldElement::<Mod631>::from_value(36_i64),
            FiniteFieldElement::<Mod631>::from_value(60_i64),
        );
        let q = EllipticCurvePoint::<FiniteFieldElement<Mod631>, CurveA30B34>::new(
            FiniteFieldElement::<Mod631>::from_value(121_i64),
            FiniteFieldElement::<Mod631>::from_value(387_i64),
        );
        let s = EllipticCurvePoint::<FiniteFieldElement<Mod631>, CurveA30B34>::new(
            FiniteFieldElement::<Mod631>::from_value(0_i64),
            FiniteFieldElement::<Mod631>::from_value(36_i64),
        );
        let order = 5.to_bigint().unwrap();

        let (fp_qs, _) = miller(&p, &(q.add_ref(&s)), &order);
        let (fp_s, _) = miller(&p, &s, &order);
        assert_eq!(fp_qs.sanitize().value, 103_i32.to_bigint().unwrap());
        assert_eq!(fp_s.sanitize().value, 219_i32.to_bigint().unwrap());
        assert_eq!(
            (fp_qs / fp_s).sanitize().value,
            473_i32.to_bigint().unwrap()
        );

        let (fq_ps, _) = miller(&q, &p.add_ref(&s.clone().neg()), &order);
        let (fq_s, _) = miller(&q, &(-s.clone()), &order);
        assert_eq!(fq_ps.sanitize().value, 284_i32.to_bigint().unwrap());
        assert_eq!(fq_s.sanitize().value, 204_i32.to_bigint().unwrap());
        assert_eq!((fq_ps / fq_s).sanitize().value, 88_i32.to_bigint().unwrap());

        let w = weil_pairing(&p, &q, &order, Some(&s));
        assert_eq!(w.sanitize().value, 242.to_bigint().unwrap());

        let p_prime = EllipticCurvePoint::<FiniteFieldElement<Mod631>, CurveA30B34>::new(
            FiniteFieldElement::<Mod631>::from_value(617_i64),
            FiniteFieldElement::<Mod631>::from_value(5_i64),
        );
        let q_prime = EllipticCurvePoint::<FiniteFieldElement<Mod631>, CurveA30B34>::new(
            FiniteFieldElement::<Mod631>::from_value(121_i64),
            FiniteFieldElement::<Mod631>::from_value(244_i64),
        );

        let w_prime = weil_pairing(&p_prime, &q_prime, &order, Some(&s));
        assert_eq!(w_prime.sanitize().value, 512_i32.to_bigint().unwrap());

        let y_prime = weil_pairing(&p_prime, &p_prime, &order, Some(&s));
        assert_eq!(y_prime.sanitize().value, 1_i32.to_bigint().unwrap());

        let y_prime2 = weil_pairing(
            &(p_prime.clone() * 3_i32.to_bigint().unwrap()),
            &(p_prime.clone() * 2_i32.to_bigint().unwrap()),
            &order,
            Some(&s),
        );
        assert_eq!(y_prime2.sanitize().value, 1_i32.to_bigint().unwrap());

        let y_prime2 = weil_pairing(
            &(p_prime.clone() * 3647912234556789907_i64.to_bigint().unwrap()),
            &(p_prime.clone() * 289_i32.to_bigint().unwrap()),
            &order,
            Some(&s),
        );
        assert_eq!(y_prime2.sanitize().value, 1_i32.to_bigint().unwrap());

        assert_eq!(p.clone() * 3_i32.to_bigint().unwrap(), p_prime.clone());
        assert_eq!(q.clone() * 4_i32.to_bigint().unwrap(), q_prime.clone());
        assert_eq!(w.pow(12_i32).sanitize(), w_prime.sanitize());
    }

    #[test]
    fn test_general_tate_pairing() {
        let p = EllipticCurvePoint::<FiniteFieldElement<Mod631>, CurveA30B34>::new(
            FiniteFieldElement::<Mod631>::from_value(36_i64),
            FiniteFieldElement::<Mod631>::from_value(60_i64),
        );
        let q = EllipticCurvePoint::<FiniteFieldElement<Mod631>, CurveA30B34>::new(
            FiniteFieldElement::<Mod631>::from_value(121_i64),
            FiniteFieldElement::<Mod631>::from_value(387_i64),
        );
        let p_prime = EllipticCurvePoint::<FiniteFieldElement<Mod631>, CurveA30B34>::new(
            FiniteFieldElement::<Mod631>::from_value(617_i64),
            FiniteFieldElement::<Mod631>::from_value(5_i64),
        );
        let q_prime = EllipticCurvePoint::<FiniteFieldElement<Mod631>, CurveA30B34>::new(
            FiniteFieldElement::<Mod631>::from_value(121_i64),
            FiniteFieldElement::<Mod631>::from_value(244_i64),
        );
        let s = EllipticCurvePoint::<FiniteFieldElement<Mod631>, CurveA30B34>::new(
            FiniteFieldElement::<Mod631>::from_value(0_i64),
            FiniteFieldElement::<Mod631>::from_value(36_i64),
        );
        let order = 5.to_bigint().unwrap();

        let tate = general_tate_pairing(&p, &q, &order, &BigInt::from(631), Some(&s));

        let tate_prime =
            general_tate_pairing(&p_prime, &q_prime, &order, &BigInt::from(631), Some(&s));

        assert_eq!(tate.pow(12_i32).sanitize(), tate_prime.sanitize());
    }

    #[test]
    fn test_tate_pairing() {
        let p = EllipticCurvePoint::<FiniteFieldElement<Mod631>, CurveA30B34>::new(
            FiniteFieldElement::<Mod631>::from_value(36_i64),
            FiniteFieldElement::<Mod631>::from_value(60_i64),
        );
        let q = EllipticCurvePoint::<FiniteFieldElement<Mod631>, CurveA30B34>::new(
            FiniteFieldElement::<Mod631>::from_value(121_i64),
            FiniteFieldElement::<Mod631>::from_value(387_i64),
        );
        let p_prime = EllipticCurvePoint::<FiniteFieldElement<Mod631>, CurveA30B34>::new(
            FiniteFieldElement::<Mod631>::from_value(617_i64),
            FiniteFieldElement::<Mod631>::from_value(5_i64),
        );
        let q_prime = EllipticCurvePoint::<FiniteFieldElement<Mod631>, CurveA30B34>::new(
            FiniteFieldElement::<Mod631>::from_value(121_i64),
            FiniteFieldElement::<Mod631>::from_value(244_i64),
        );
        let order = 5.to_bigint().unwrap();

        let tate = tate_pairing(&p, &q, &order, &BigInt::from(631));

        let tate_prime = tate_pairing(&p_prime, &q_prime, &order, &BigInt::from(631));

        assert_eq!(tate.pow(12_i32).sanitize(), tate_prime.sanitize());
    }
}
