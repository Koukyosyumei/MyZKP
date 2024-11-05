use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::fmt::Debug;
use std::marker::PhantomData;
use std::ops::{Add, Mul, Neg, Sub};

use crate::modules::field::Field;

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
    fn new(x: F, y: F) -> Self {
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

    pub fn line_slope(&self, other: Self) -> F {
        let a = F::from_value(E::get_a());
        // let b = E::get_b();

        let x1 = self.x.clone().unwrap();
        let y1 = self.y.clone().unwrap();
        let x2 = other.x.clone().unwrap();
        let y2 = other.y.clone().unwrap();

        if self.x.clone() == other.x.clone() {
            ((x1.clone() * x1.clone()) * (3_i64) + a.clone()) / (y1.clone() * (2_i64))
        } else {
            (y2.clone() - y1.clone()) / (x2.clone() - x1.clone())
        }
    }
}

impl<F: Field, E: EllipticCurve> Add for EllipticCurvePoint<F, E> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        if self.is_point_at_infinity() {
            return other;
        }
        if other.is_point_at_infinity() {
            return self;
        }

        let m = self.line_slope(other.clone());

        // If x coordinates are the same
        if self.x == other.x {
            // Check if they are inverses of each other
            if self.y != other.y {
                // TODO: check self.y = -other.y
                return EllipticCurvePoint::point_at_infinity();
            } else {
                // P = Q, handle point doubling
                let x1 = self.x.clone().unwrap();
                let y1 = self.y.clone().unwrap();

                // let m = ((x1.clone() * x1.clone()).mul_scalar(3_i64) + self.curve.a.clone())
                //    / (y1.clone().mul_scalar(2_i64));

                let x3 = m.clone() * m.clone() - x1.clone() - x1.clone();
                let y3 = m * (x1 - x3.clone()) - y1;

                return EllipticCurvePoint::new(x3, y3);
            }
        } else {
            // P != Q, handle regular addition
            let x1 = self.x.clone().unwrap();
            let y1 = self.y.clone().unwrap();
            let x2 = other.x.clone().unwrap();
            // let y2 = other.y.clone().unwrap();

            // let m = (y2.clone() - y1.clone()) / (x2.clone() - x1.clone());
            let x3 = m.clone() * m.clone() - x1.clone() - x2.clone();
            let y3 = m * (x1 - x3.clone()) - y1;

            return EllipticCurvePoint::new(x3, y3);
        }
    }
}

impl<F: Field, E: EllipticCurve> Mul<BigInt> for EllipticCurvePoint<F, E> {
    type Output = Self;

    fn mul(self, scalar: BigInt) -> Self {
        if scalar.is_zero() {
            // Return the point at infinity for scalar * 0
            return EllipticCurvePoint::point_at_infinity();
        }

        let mut result = EllipticCurvePoint::point_at_infinity();
        let mut current = self.clone(); // Start with the current point
        let mut scalar_bits = scalar.clone();

        while scalar_bits > BigInt::zero() {
            if &scalar_bits & BigInt::one() == BigInt::one() {
                result = result + current.clone(); // Add the current point if the bit is 1
            }
            current = current.clone() + current.clone(); // Double the point
            scalar_bits >>= 1; // Move to the next bit
        }

        result
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

pub fn get_lambda<F: Field, E: EllipticCurve>(
    p: EllipticCurvePoint<F, E>,
    q: EllipticCurvePoint<F, E>,
    r: EllipticCurvePoint<F, E>,
) -> F {
    let p_x = p.x.clone().unwrap();
    let p_y = p.y.clone().unwrap();
    let q_x = q.x.clone().unwrap();
    // let q_y = q.y.clone().unwrap();
    let r_x = r.x.clone().unwrap();
    let r_y = r.y.clone().unwrap();

    if (p == q && p_y.clone() == F::zero()) || (p != q && p_x.clone() == q_x.clone()) {
        return r_x.clone() - p_x.clone();
    }
    let slope = p.line_slope(q.clone());
    let numerator = r_y.clone() - p_y.clone() - slope.clone() * (r_x.clone() - p_x.clone());
    let denominator = r_x.clone() + p_x.clone() + q_x.clone() - slope.clone() * slope.clone();
    return numerator / denominator;
}

pub fn miller<F: Field, E: EllipticCurve>(
    p: EllipticCurvePoint<F, E>,
    q: EllipticCurvePoint<F, E>,
    m: BigInt,
) -> F {
    if p == q {
        F::one();
    }

    let mut f = F::one();
    let mut t = p.clone();

    for i in (1..m.bits()).rev() {
        f = (f.clone() * f.clone()) * (get_lambda(t.clone(), t.clone(), q.clone()));
        t = t.clone() + t.clone();
        if m.bit(i) {
            f = f * (get_lambda(t.clone(), p.clone(), q.clone()));
            t = t.clone() + p.clone();
        }
    }

    f
}

pub fn weil_pairing<F: Field, E: EllipticCurve>(
    p: EllipticCurvePoint<F, E>,
    q: EllipticCurvePoint<F, E>,
    m: BigInt,
    s: Option<EllipticCurvePoint<F, E>>,
) -> F {
    let s_value = s.unwrap();
    let fp_qs = miller(p.clone(), q.clone() + s_value.clone(), m.clone());
    let fp_s = miller(p.clone(), s_value.clone(), m.clone());
    let fq_ps = miller(q.clone(), p.clone() - s_value.clone(), m.clone());
    let fq_s = miller(q.clone(), -s_value.clone(), m.clone());

    return (fp_qs / fp_s) / (fq_ps / fq_s);
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
    use crate::{
        define_myzkp_modulus_type,
        modules::field::{FiniteFieldElement, ModulusValue},
    };
    use num_bigint::{BigInt, ToBigInt};
    use std::str::FromStr;

    #[test]
    fn test_weil_pairing() {
        define_myzkp_modulus_type!(Mod631, "631");
        define_myzkp_curve_type!(CurveA30B34, "30", "34");

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

        let fp_qs = miller(p.clone(), q.clone() + s.clone(), order.clone());
        let fp_s = miller(p.clone(), s.clone(), order.clone());
        assert_eq!(fp_qs.sanitize().value, 103_i32.to_bigint().unwrap());
        assert_eq!(fp_s.sanitize().value, 219_i32.to_bigint().unwrap());
        assert_eq!(
            (fp_qs / fp_s).sanitize().value,
            473_i32.to_bigint().unwrap()
        );

        let fq_ps = miller(q.clone(), p.clone() - s.clone(), order.clone());
        let fq_s = miller(q.clone(), -s.clone(), order.clone());
        assert_eq!(fq_ps.sanitize().value, 284_i32.to_bigint().unwrap());
        assert_eq!(fq_s.sanitize().value, 204_i32.to_bigint().unwrap());
        assert_eq!((fq_ps / fq_s).sanitize().value, 88_i32.to_bigint().unwrap());

        let w = weil_pairing(p.clone(), q.clone(), order.clone(), Some(s.clone()));
        assert_eq!(w.sanitize().value, 242.to_bigint().unwrap());

        let p_prime = EllipticCurvePoint::<FiniteFieldElement<Mod631>, CurveA30B34>::new(
            FiniteFieldElement::<Mod631>::from_value(617_i64),
            FiniteFieldElement::<Mod631>::from_value(5_i64),
        );
        let q_prime = EllipticCurvePoint::<FiniteFieldElement<Mod631>, CurveA30B34>::new(
            FiniteFieldElement::<Mod631>::from_value(121_i64),
            FiniteFieldElement::<Mod631>::from_value(244_i64),
        );

        let w_prime = weil_pairing(
            p_prime.clone(),
            q_prime.clone(),
            order.clone(),
            Some(s.clone()),
        );
        assert_eq!(w_prime.sanitize().value, 512_i32.to_bigint().unwrap());

        assert_eq!(p.clone() * 3_i32.to_bigint().unwrap(), p_prime.clone());
        assert_eq!(q.clone() * 4_i32.to_bigint().unwrap(), q_prime.clone());
        assert_eq!(
            w.pow(12_i32.to_bigint().unwrap()).sanitize(),
            w_prime.sanitize()
        );
    }
}
