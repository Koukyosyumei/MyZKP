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
