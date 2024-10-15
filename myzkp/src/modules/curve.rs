use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_traits::{One, Signed, Zero};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct EllipticCurvePoint {
    x: BigInt,       // x-coordinate of the elliptic curve point
    y: BigInt,       // y-coordinate of the elliptic curve point
    modulus: BigInt, // The modulus (the prime field over which the curve is defined)
    a: BigInt,       // The 'a' parameter of the elliptic curve equation y^2 = x^3 + ax + b
    b: BigInt,       // The 'b' parameter of the elliptic curve equation y^2 = x^3 + ax + b
}

const DEFAULT_K_MODULES: i64 = 3_i64 * (1 << 30) + 1;
