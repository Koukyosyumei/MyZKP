use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_traits::{One, Signed, Zero};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::modules::field::Field;
use crate::modules::field::FiniteFieldElement;
use crate::modules::polynomial::Polynomial;

#[derive(Debug, Clone, PartialEq)]
pub struct EllipticCurve<F: Field> {
    pub a: F,
    pub b: F,
}

#[derive(Debug, Clone, PartialEq)]
pub struct EllipticCurvePoint<F: Field> {
    pub x: Option<F>,
    pub y: Option<F>,
    pub curve: EllipticCurve<F>,
}

impl<F: Field> EllipticCurvePoint<F> {
    fn new(x: F, y: F, curve: EllipticCurve<F>) -> Self {
        EllipticCurvePoint {
            x: Some(x),
            y: Some(y),
            curve: curve,
        }
    }

    pub fn point_at_infinity(curve: EllipticCurve<F>) -> Self {
        EllipticCurvePoint {
            x: None,
            y: None,
            curve: curve,
        }
    }

    pub fn is_point_at_infinity(&self) -> bool {
        self.x.is_none() || self.y.is_none()
    }

    pub fn inverse(&self) -> Self {
        if self.is_point_at_infinity() {
            return self.clone();
        }

        EllipticCurvePoint::new(
            self.x.clone().unwrap(),
            -(self.y.clone().unwrap()),
            self.curve.clone(),
        )
    }

    pub fn line_slope(&self, other: Self) -> F {
        let x1 = self.x.clone().unwrap();
        let y1 = self.y.clone().unwrap();
        let x2 = other.x.clone().unwrap();
        let y2 = other.y.clone().unwrap();

        if self.x.clone() == other.x.clone() {
            ((x1.clone() * x1.clone()).mul_scalar(3_i64) + self.curve.a.clone())
                / (y1.clone().mul_scalar(2_i64))
        } else {
            (y2.clone() - y1.clone()) / (x2.clone() - x1.clone())
        }
    }

    pub fn mul_scalar(&self, scalar: BigInt) -> Self {
        if scalar.is_zero() {
            // Return the point at infinity for scalar * 0
            return EllipticCurvePoint::point_at_infinity(self.curve.clone());
        }

        let mut result = EllipticCurvePoint::point_at_infinity(self.curve.clone());
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

impl<F: Field> Add for EllipticCurvePoint<F> {
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
                return EllipticCurvePoint::point_at_infinity(self.curve.clone());
            } else {
                // P = Q, handle point doubling
                let x1 = self.x.clone().unwrap();
                let y1 = self.y.clone().unwrap();

                // let m = ((x1.clone() * x1.clone()).mul_scalar(3_i64) + self.curve.a.clone())
                //    / (y1.clone().mul_scalar(2_i64));

                let x3 = m.clone() * m.clone() - x1.clone() - x1.clone();
                let y3 = m * (x1 - x3.clone()) - y1;

                return EllipticCurvePoint::new(x3, y3, self.curve.clone());
            }
        } else {
            // P != Q, handle regular addition
            let x1 = self.x.clone().unwrap();
            let y1 = self.y.clone().unwrap();
            let x2 = other.x.clone().unwrap();
            let y2 = other.y.clone().unwrap();

            // let m = (y2.clone() - y1.clone()) / (x2.clone() - x1.clone());
            let x3 = m.clone() * m.clone() - x1.clone() - x2.clone();
            let y3 = m * (x1 - x3.clone()) - y1;

            return EllipticCurvePoint::new(x3, y3, self.curve.clone());
        }
    }
}

impl<F: Field> Neg for EllipticCurvePoint<F> {
    type Output = Self;
    fn neg(self) -> Self {
        if self.is_point_at_infinity() {
            return self;
        } else {
            EllipticCurvePoint::new(
                self.x.clone().unwrap(),
                -self.y.clone().unwrap(),
                self.curve,
            )
        }
    }
}

impl<F: Field> Sub for EllipticCurvePoint<F> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self + (-other)
    }
}

pub fn get_lambda<F: Field>(
    p: EllipticCurvePoint<F>,
    q: EllipticCurvePoint<F>,
    r: EllipticCurvePoint<F>,
    modulus: Option<BigInt>,
) -> F {
    let p_x = p.x.clone().unwrap();
    let p_y = p.y.clone().unwrap();
    let q_x = q.x.clone().unwrap();
    let q_y = q.y.clone().unwrap();
    let r_x = r.x.clone().unwrap();
    let r_y = r.y.clone().unwrap();

    if (p == q && p_y.clone() == F::zero(modulus)) || (p != q && p_x.clone() == q_x.clone()) {
        return r_x.clone() - p_x.clone();
    }
    let slope = p.line_slope(q.clone());
    let numerator = r_y.clone() - p_y.clone() - slope.clone() * (r_x.clone() - p_x.clone());
    let denominator = r_x.clone() + p_x.clone() + q_x.clone() - slope.clone() * slope.clone();
    return numerator / denominator;
}

pub fn miller<F: Field>(
    p: EllipticCurvePoint<F>,
    q: EllipticCurvePoint<F>,
    m: BigInt,
    modulus: Option<BigInt>,
) -> F {
    if p == q {
        F::one(modulus.clone());
    }

    let mut f = F::one(modulus.clone());
    let mut t = p.clone();

    for i in (1..m.bits()).rev() {
        f = (f.clone() * f.clone())
            * (get_lambda(t.clone(), t.clone(), q.clone(), modulus.clone()));
        t = t.clone() + t.clone();
        if m.bit(i) {
            f = f * (get_lambda(t.clone(), p.clone(), q.clone(), modulus.clone()));
            t = t.clone() + p.clone();
        }
    }

    f
}

pub fn weil_pairing<F: Field>(
    p: EllipticCurvePoint<F>,
    q: EllipticCurvePoint<F>,
    m: BigInt,
    s: Option<EllipticCurvePoint<F>>,
    modulus: Option<BigInt>,
) -> F {
    let s_value = s.unwrap();
    let fp_qs = miller(
        p.clone(),
        q.clone() + s_value.clone(),
        m.clone(),
        modulus.clone(),
    );
    let fp_s = miller(p.clone(), s_value.clone(), m.clone(), modulus.clone());
    let fq_qs = miller(
        q.clone(),
        q.clone() - s_value.clone(),
        m.clone(),
        modulus.clone(),
    );
    let fq_s = miller(q.clone(), -s_value.clone(), m.clone(), modulus.clone());

    return (fp_qs / fp_s) * (fq_qs / fq_s);
}

#[test]
fn test_weil_pairing() {
    let modulus = Some(631_i64.to_bigint().unwrap());
    let a = FiniteFieldElement::new(30_i64.to_bigint().unwrap(), modulus.clone());
    let b = FiniteFieldElement::new(34_i64.to_bigint().unwrap(), modulus.clone());
    let curve = EllipticCurve { a, b };

    let p = EllipticCurvePoint::new(
        FiniteFieldElement::new(36_i64.to_bigint().unwrap(), modulus.clone()),
        FiniteFieldElement::new(60_i64.to_bigint().unwrap(), modulus.clone()),
        curve.clone(),
    );
    let q = EllipticCurvePoint::new(
        FiniteFieldElement::new(121_i64.to_bigint().unwrap(), modulus.clone()),
        FiniteFieldElement::new(387_i64.to_bigint().unwrap(), modulus.clone()),
        curve.clone(),
    );
    let s = EllipticCurvePoint::new(
        FiniteFieldElement::new(0_i64.to_bigint().unwrap(), modulus.clone()),
        FiniteFieldElement::new(36_i64.to_bigint().unwrap(), modulus.clone()),
        curve.clone(),
    );
    let order = 5.to_bigint().unwrap();

    let fp_qs = miller(
        p.clone(),
        q.clone() + s.clone(),
        order.clone(),
        modulus.clone(),
    );
    let fp_s = miller(p.clone(), s.clone(), order.clone(), modulus.clone());
    assert_eq!(fp_qs.value, 103.to_bigint().unwrap());
    assert_eq!(fp_s.value, 219.to_bigint().unwrap());
    assert_eq!((fp_qs / fp_s).value, 473.to_bigint().unwrap());
}
