use std::str::FromStr;

use lazy_static::lazy_static;
use num_bigint::BigInt;
use num_bigint::ToBigInt;
use num_traits::{One, Zero};
use paste::paste;

use crate::modules::algebra::curve::curve::{
    get_lambda, miller, EllipticCurve, EllipticCurvePoint,
};
use crate::modules::algebra::efield::{ExtendedFieldElement, IrreduciblePoly};
use crate::modules::algebra::field::{FiniteFieldElement, ModEIP197, ModulusValue};
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::ring::Ring;
use crate::{define_myzkp_curve_type, define_myzkp_modulus_type};

define_myzkp_modulus_type!(
    BN128Modulus,
    "21888242871839275222246405745257275088696311157297823662689037894645226208583"
);
define_myzkp_curve_type!(BN128Curve, "0", "3");

lazy_static! {
    static ref ATE_LOOP_COUNT: BigInt = BigInt::from_str("29793968203157093288").unwrap();
}

pub type Fq = FiniteFieldElement<BN128Modulus>;
pub type FqOrder = FiniteFieldElement<ModEIP197>;
pub type G1Point = EllipticCurvePoint<Fq, BN128Curve>;

// Define Fq2 as a quadratic extension field
#[derive(Debug, Clone, PartialEq, Hash)]
pub struct Fq2Poly;
lazy_static! {
    // x^2 + 1
    static ref MODULUS_Fq2: Polynomial<Fq> = Polynomial {
        coef: vec![Fq::one(), Fq::zero(), Fq::one()]
    };
}

impl IrreduciblePoly<Fq> for Fq2Poly {
    fn modulus() -> &'static Polynomial<Fq> {
        &MODULUS_Fq2
    }
}
pub type Fq2 = ExtendedFieldElement<BN128Modulus, Fq2Poly>;
pub type G2Point = EllipticCurvePoint<Fq2, BN128Curve>;

// Define Fq2 as a quadratic extension field
#[derive(Debug, Clone, PartialEq, Hash)]
pub struct Fq12Poly;
lazy_static! {
    // x^2 + 1
    static ref MODULUS_Fq12: Polynomial<Fq> = Polynomial {
        // x^12 -18x^6 + 82
        coef: vec![
            Fq::from_value(82_i32),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::from_value(-18_i32),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::one(),
        ],
    };
}
impl IrreduciblePoly<Fq> for Fq12Poly {
    fn modulus() -> &'static Polynomial<Fq> {
        &MODULUS_Fq12
    }
}
type Fq12 = ExtendedFieldElement<BN128Modulus, Fq12Poly>;
type G12Point = EllipticCurvePoint<Fq12, BN128Curve>;

pub fn cast_g1_to_g12(g: G1Point) -> G12Point {
    if g.is_point_at_infinity() {
        return G12Point::point_at_infinity();
    }

    G12Point::new(
        Fq12::new(Polynomial {
            coef: vec![g.x.unwrap()],
        }),
        Fq12::new(Polynomial {
            coef: vec![g.y.unwrap()],
        }),
    )
}

pub fn twist_g2_to_g12(g: G2Point) -> G12Point {
    if g.is_point_at_infinity() {
        return G12Point::point_at_infinity();
    }

    let w = Fq12::new(Polynomial {
        coef: vec![Fq::zero(), Fq::one()],
    });

    let x = g.x.unwrap();
    let y = g.y.unwrap();
    let x_coeff = vec![
        x.poly.coef[0].clone() - x.poly.coef[1].clone() * Fq::from_value(9_i32),
        x.poly.coef[1].clone(),
    ];
    let y_coeff = vec![
        y.poly.coef[0].clone() - y.poly.coef[1].clone() * Fq::from_value(9_i32),
        y.poly.coef[1].clone(),
    ];

    let nx = Fq12::new(Polynomial {
        coef: vec![
            x_coeff[0].clone(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            x_coeff[1].clone(),
        ],
    });
    let ny = Fq12::new(Polynomial {
        coef: vec![
            y_coeff[0].clone(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            y_coeff[1].clone(),
        ],
    });

    return G12Point::new(
        nx * (w.pow(2_i32.to_bigint().unwrap())),
        ny * (w.pow(3_i32.to_bigint().unwrap())),
    );
}

pub fn optimal_ate_pairing(p_g1: &G1Point, q_g2: &G2Point) -> Fq12 {
    // https://eprint.iacr.org/2010/354.pdf
    let p: G12Point = cast_g1_to_g12(p_g1.clone());
    let q: G12Point = twist_g2_to_g12(q_g2.clone());
    let m = BN128Modulus::modulus();

    if p.is_point_at_infinity() || q.is_point_at_infinity() {
        return Fq12::one();
    }

    let mut f = Fq12::one();
    if p != q {
        let out = miller(&q, &p, &ATE_LOOP_COUNT);
        f = out.0;
        let mut r = out.1;
        // Assert: r == multiply(&q, &ate_loop_count)

        let q1 = G12Point::new(
            q.x.clone().unwrap().pow(m.clone()),
            q.y.clone().unwrap().pow(m.clone()),
        );

        let nq2 = G12Point::new(
            q1.x.clone().unwrap().pow(m.clone()),
            -q1.y.clone().unwrap().pow(m.clone()),
        );

        f *= get_lambda(&r, &q1, &p);
        r += q1;
        f *= get_lambda(&r, &nq2, &p);
    }

    let exp = (m.pow(12) - BigInt::one()) / (BN128::order());
    f.pow(exp)
}

pub struct BN128;

impl BN128 {
    pub fn generator_g1() -> G1Point {
        G1Point::new(Fq::from_value(1), Fq::from_value(2))
    }

    pub fn generator_g2() -> G2Point {
        G2Point::new(
            Fq2::new(Polynomial {
                coef: vec![
                    Fq::from_value(BigInt::from_str("10857046999023057135944570762232829481370756359578518086990519993285655852781").unwrap()),
                    Fq::from_value(BigInt::from_str("11559732032986387107991004021392285783925812861821192530917403151452391805634").unwrap()),
                ],
            }),
            Fq2::new(Polynomial {
                coef: vec![
                    Fq::from_value(BigInt::from_str("8495653923123431417604973247489272438418190587263600148770280649306958101930").unwrap()),
                    Fq::from_value(BigInt::from_str("4082367875863433681332203403145435568316851327593401208105741076214120093531").unwrap()),
                ],
            }),
        )
    }

    pub fn order() -> BigInt {
        BigInt::from_str(
            "21888242871839275222246405745257275088548364400416034343698204186575808495617",
        )
        .unwrap()
    }

    pub fn get_b() -> Fq {
        Fq::from_value(3_i32)
    }

    pub fn get_b2() -> Fq2 {
        Fq2::new(Polynomial {
            coef: vec![Fq::from_value(3_i32)],
        }) / Fq2::new(Polynomial {
            coef: vec![Fq::from_value(9_i32), Fq::from_value(1_i32)],
        })
    }

    pub fn get_b12() -> Fq12 {
        Fq12::new(Polynomial {
            coef: vec![Fq::from_value(3_i32)],
        })
    }
}

// Test the implementation
#[cfg(test)]
mod tests {
    use super::*;

    use crate::modules::algebra::field::Field;

    #[test]
    fn test_fq() {
        assert_eq!(Fq::from_value(2) * Fq::from_value(2), Fq::from_value(4));
        assert_eq!(
            Fq::from_value(2) / Fq::from_value(7) + Fq::from_value(9) / Fq::from_value(7),
            Fq::from_value(11) / Fq::from_value(7)
        );
        assert_eq!(
            Fq::from_value(2) * Fq::from_value(7) + Fq::from_value(9) * Fq::from_value(7),
            Fq::from_value(11) * Fq::from_value(7)
        );
    }

    #[test]
    fn test_fq2() {
        let x = Fq2::new(Polynomial {
            coef: vec![Fq::one(), Fq::zero()],
        });
        let f = Fq2::new(Polynomial {
            coef: vec![Fq::one(), Fq::from_value(2_i32)],
        });
        let fpx = Fq2::new(Polynomial {
            coef: vec![Fq::from_value(2_i32), Fq::from_value(2_i32)],
        });
        let fpx2 = Fq2::new(Polynomial {
            coef: vec![Fq::from_value(2_i32), Fq::from_value(1_i32)],
        });
        let one = Fq2::one();
        assert_eq!(x.add_ref(&f), fpx);
        assert_eq!(fpx2.div_ref(&fpx2), one);
        assert_eq!(
            one.div_ref(&f) + x.div_ref(&f),
            (one.clone() + x.clone()) / f.clone()
        );
        assert_eq!(
            one.clone() * f.clone() + x.clone() * f.clone(),
            (one.clone() + x.clone()) * f.clone()
        );
        assert_eq!(
            x.clone()
                .pow(BN128Modulus::modulus().pow(2) - BigInt::one()),
            one
        );
    }

    #[test]
    fn test_g1() {
        let g1 = BN128::generator_g1();
        assert_eq!(
            g1.clone().y.unwrap().pow(2) - g1.clone().x.unwrap().pow(3),
            BN128::get_b()
        );
        assert_eq!(
            g1.clone() * 2 + g1.clone() + g1.clone(),
            (g1.clone() * 2) * 2
        );
        assert_eq!(
            g1.clone() * 9 + g1.clone() * 5,
            g1.clone() * 12 + g1.clone() * 2,
        );
        assert!((g1.clone() * BN128::order()).is_point_at_infinity());
    }

    #[test]
    fn test_g2() {
        let g2 = BN128::generator_g2();
        assert_eq!(
            g2.clone().y.unwrap().pow(2) - g2.clone().x.unwrap().pow(3),
            BN128::get_b2()
        );
        assert_eq!(
            g2.clone() * 2 + g2.clone() + g2.clone(),
            (g2.clone() * 2) * 2
        );
        assert_eq!(
            g2.clone() * 9 + g2.clone() * 5,
            g2.clone() * 12 + g2.clone() * 2,
        );
        assert!((g2.clone() * BN128::order()).is_point_at_infinity());
    }

    #[test]
    fn test_g12() {
        let g2 = BN128::generator_g2();
        let g12 = twist_g2_to_g12(g2);
        let g12_9 = g12.clone() * 9;
        assert_eq!(
            g12_9.clone().y.unwrap().pow(2) - g12_9.clone().x.unwrap().pow(3),
            BN128::get_b12()
        );
        assert_eq!(
            g12.clone() * 2 + g12.clone() + g12.clone(),
            (g12.clone() * 2) * 2
        );
        assert_eq!(
            g12.clone() * 9 + g12.clone() * 5,
            g12.clone() * 12 + g12.clone() * 2,
        );
        assert!((g12.clone() * BN128::order()).is_point_at_infinity());
    }

    #[test]
    fn test_pairing() {
        let g1 = BN128::generator_g1();
        let g2 = BN128::generator_g2();
        let p1 = optimal_ate_pairing(&g1.clone(), &g2.clone());
        let pn1 = optimal_ate_pairing(&-g1.clone(), &g2.clone());
        assert_eq!(p1.clone() * pn1.clone(), Fq12::one());
        let np1 = optimal_ate_pairing(&g1.clone(), &-g2.clone());
        assert_eq!(p1.clone() * np1.clone(), Fq12::one());
        assert_eq!(pn1.clone(), np1.clone());
        let p2 = optimal_ate_pairing(&(g1.mul_ref(2)), &g2.clone());
        assert_eq!(p1.clone() * p1.clone(), p2);
        assert!(p1 != p2);
        assert!(p1 != np1);
        assert!(p2 != np1);
        let po2 = optimal_ate_pairing(&g1.clone(), &(g2.mul_ref(2)));
        assert_eq!(p1.clone() * p1.clone(), po2);
        let p3 = optimal_ate_pairing(&(g1.mul_ref(37)), &(g2.mul_ref(27)));
        let po3 = optimal_ate_pairing(&(g1 * 999), &g2);
        assert_eq!(p3, po3);
    }
}
