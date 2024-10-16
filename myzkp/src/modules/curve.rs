use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_traits::{One, Signed, Zero};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::modules::field::Field;
use crate::modules::field::FiniteFieldElement;

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

                let m = ((x1.clone() * x1.clone()).mul_scalar(3_i64) + self.curve.a.clone())
                    / (y1.clone().mul_scalar(2_i64));

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

            let m = (y2.clone() - y1.clone()) / (x2.clone() - x1.clone());
            let x3 = m.clone() * m.clone() - x1.clone() - x2.clone();
            let y3 = m * (x1 - x3.clone()) - y1;

            return EllipticCurvePoint::new(x3, y3, self.curve.clone());
        }
    }
}

/*
impl<F: Field> Sub for EllipticCurvePoint<F> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self.add(&other.neg())
    }
}

impl<F: Field> Neg for EllipticCurvePoint<F> {
    type Output = Self;

    fn neg(self) -> Self {
        self.neg()
    }
}

impl<F: Field> Mul<BigInt> for EllipticCurvePoint<F> {
    type Output = Self;

    fn mul(self, scalar: BigInt) -> Self {
        self.mul_scalar(scalar)
    }
}
    */

/*
use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_traits::{One, Signed, Zero};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::modules::field::Field;
use crate::modules::field::FiniteFieldElement;
use crate::modules::polynomial::Polynomial;

impl<F: Field> EllipticCurvePoint<F> {
    // ... (existing methods)

    pub fn scalar_mul(&self, scalar: &BigInt) -> Self {
        let mut result = EllipticCurvePoint::point_at_infinity(self.curve.clone());
        let mut temp = self.clone();
        let mut n = scalar.clone();

        while n > BigInt::zero() {
            if n.is_odd() {
                result = result + temp.clone();
            }
            temp = temp.clone() + temp.clone();
            n >>= 1;
        }

        result
    }

    pub fn miller_function(&self, q: &Self, m: &BigInt) -> Polynomial<F> {
        let mut f = Polynomial::from_constant(F::one(None));
        let mut t = self.clone();

        for i in (1..m.bits()).rev() {
            let line = line_function(&t, &t);
            let vertical = vertical_line(&(t.clone() + t.clone()));
            f = f.clone() * f.clone() * line / vertical;

            if m.bit(i) {
                let line = line_function(&t, self);
                let vertical = vertical_line(&(t.clone() + self.clone()));
                f = f * line / vertical;
                t = t + self.clone();
            }

            t = t + t;
        }

        f.eval(&q.x.unwrap()) / f.eval(&(q + self.inverse()).x.unwrap())
    }
}

fn line_function<F: Field>(p: &EllipticCurvePoint<F>, q: &EllipticCurvePoint<F>) -> Polynomial<F> {
    if p.is_point_at_infinity() || q.is_point_at_infinity() {
        return Polynomial::from_constant(F::one(None));
    }

    if p == q {
        let m = ((p.x.clone().unwrap() * p.x.clone().unwrap()).mul_scalar(3_i64) + p.curve.a.clone())
            / (p.y.clone().unwrap().mul_scalar(2_i64));
        let c = p.y.clone().unwrap() - m.clone() * p.x.clone().unwrap();
        Polynomial::from_coefficients(vec![c, m, F::one(None).neg()])
    } else {
        let m = (q.y.clone().unwrap() - p.y.clone().unwrap()) / (q.x.clone().unwrap() - p.x.clone().unwrap());
        let c = p.y.clone().unwrap() - m.clone() * p.x.clone().unwrap();
        Polynomial::from_coefficients(vec![c, m, F::one(None).neg()])
    }
}

fn vertical_line<F: Field>(p: &EllipticCurvePoint<F>) -> Polynomial<F> {
    if p.is_point_at_infinity() {
        Polynomial::from_constant(F::one(None))
    } else {
        Polynomial::from_coefficients(vec![p.x.clone().unwrap().neg(), F::one(None)])
    }
}

pub fn weil_pairing<F: Field>(p: &EllipticCurvePoint<F>, q: &EllipticCurvePoint<F>, m: &BigInt) -> F {
    let fp = p.miller_function(q, m);
    let fq = q.miller_function(p, m);
    fp / fq
}
*/
