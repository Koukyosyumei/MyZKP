use std::fmt;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Rem, Sub, SubAssign};

use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};

use crate::modules::algebra::curve::curve::{EllipticCurve, EllipticCurvePoint};
use crate::modules::algebra::field::Field;

pub fn tensor_product<F: Field>(a: &Vec<F>, b: &Vec<F>) -> Vec<F> {
    // a * b^T
    let ab = a
        .iter()
        .map(|v| b.iter().map(|u| v.mul_ref(u)).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    (0..b.len())
        .map(|i| ab.iter().map(|vs| vs[i].clone()).collect::<Vec<_>>())
        .flatten()
        .collect::<Vec<_>>()
}

mod tests {
    use super::*;

    use std::str::FromStr;

    use lazy_static::lazy_static;
    use num_bigint::BigInt;
    use num_bigint::ToBigInt;
    use paste::paste;
    use serde::{Deserialize, Serialize};

    use crate::define_myzkp_modulus_type;
    use crate::modules::algebra::field::{FiniteFieldElement, ModulusValue};
    use crate::modules::algebra::polynomial::Polynomial;
    use crate::modules::algebra::ring::Ring;

    define_myzkp_modulus_type!(Mod31, "31");
    type F31 = FiniteFieldElement<Mod31>;

    #[test]
    fn test_tensor_product_1() {
        let a = vec![F31::from_value(2), F31::from_value(3)];
        let b = vec![F31::from_value(4), F31::from_value(5), F31::from_value(6)];
        let c = tensor_product(&a, &b);
        let c_groundtruth = vec![
            F31::from_value(8),
            F31::from_value(12),
            F31::from_value(10),
            F31::from_value(15),
            F31::from_value(12),
            F31::from_value(18),
        ];
        assert_eq!(c.len(), c_groundtruth.len());
        for (u, v) in c.iter().zip(c_groundtruth.iter()) {
            assert_eq!(u, v);
        }
    }

    #[test]
    fn test_tensor_product_2() {
        let x1 = vec![F31::from_value(1), F31::from_value(2)];
        let x2 = vec![F31::from_value(1), F31::from_value(3)];
        let x3 = vec![F31::from_value(1), F31::from_value(4)];
        let c = tensor_product(&tensor_product(&x1, &x2), &x3);
        let c_groundtruth = vec![
            F31::from_value(1),
            F31::from_value(2),
            F31::from_value(3),
            F31::from_value(6),
            F31::from_value(4),
            F31::from_value(8),
            F31::from_value(12),
            F31::from_value(24),
        ];
        assert_eq!(c.len(), c_groundtruth.len());
        for (u, v) in c.iter().zip(c_groundtruth.iter()) {
            assert_eq!(u, v);
        }
    }
}
