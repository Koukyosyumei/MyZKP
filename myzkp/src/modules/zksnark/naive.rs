use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_traits::One;
use num_traits::Zero;
use std::str::FromStr;

use crate::modules::bn128::{optimal_ate_pairing, Fq, Fq2, FqOrder, G1Point, G2Point};
use crate::modules::curve::{EllipticCurve, EllipticCurvePoint};
use crate::modules::field::Field;
use crate::modules::polynomial::Polynomial;
use crate::modules::qap::QAP;
use crate::modules::ring::Ring;

#[derive(Debug, Clone)]
pub struct ProofKey {
    g1_ell_i_vec: Vec<G1Point>,
    g1_r_i_vec: Vec<G1Point>,
    g2_r_i_vec: Vec<G2Point>,
    g1_o_i_vec: Vec<G1Point>,
    g1_alpha_ell_i_vec: Vec<G1Point>,
    g1_alpha_r_i_vec: Vec<G1Point>,
    g1_alpha_o_i_vec: Vec<G1Point>,
    g1_sj_vec: Vec<G1Point>,
}

#[derive(Debug, Clone)]
pub struct VerificationKey {
    g2_alpha: G2Point,
    g2_t_s: G2Point,
}

#[derive(Debug, Clone)]
pub struct Proof {
    g1_ell: G1Point,
    g1_r: G1Point,
    g2_r: G2Point,
    g1_o: G1Point,
    g1_ell_prime: G1Point,
    g1_r_prime: G1Point,
    g1_o_prime: G1Point,
    g1_h: G1Point,
}

pub fn generate_challenge_vec<F1: Field, F2: Field, E: EllipticCurve>(
    point: &EllipticCurvePoint<F1, E>,
    poly_vec: &[Polynomial<F2>],
    s: &F2,
) -> Vec<EllipticCurvePoint<F1, E>> {
    poly_vec
        .iter()
        .map(|poly| point.mul_ref(poly.eval(s).sanitize().get_value()))
        .collect()
}

pub fn generate_alpha_challenge_vec<F1: Field, F2: Field, E: EllipticCurve>(
    point: &EllipticCurvePoint<F1, E>,
    poly_vec: &[Polynomial<F2>],
    s: &F2,
    alpha: &F2,
) -> Vec<EllipticCurvePoint<F1, E>> {
    poly_vec
        .iter()
        .map(|poly| point.mul_ref((alpha.mul_ref(&poly.eval(s).sanitize())).get_value()))
        .collect()
}

pub fn generate_s_powers<F1: Field, F2: Field, E: EllipticCurve>(
    point: &EllipticCurvePoint<F1, E>,
    s: &F2,
    m: usize,
) -> Vec<EllipticCurvePoint<F1, E>> {
    let mut powers = Vec::with_capacity(m + 1);
    let mut current = F2::one();
    for _ in 0..=m {
        powers.push(point.mul_ref(current.get_value()));
        current = current * s.clone();
    }
    powers
}

fn accumulate_curve_points<F1: Field, F2: Field, E: EllipticCurve>(
    g_vec: &[EllipticCurvePoint<F1, E>],
    assignment: &[F2],
) -> EllipticCurvePoint<F1, E> {
    g_vec.iter().zip(assignment.iter()).fold(
        EllipticCurvePoint::<F1, E>::point_at_infinity(),
        |acc, (g, &ref a)| acc + g.mul_ref(a.get_value()),
    )
}

fn accumulate_polynomials<F: Field>(poly_vec: &[Polynomial<F>], assignment: &[F]) -> Polynomial<F> {
    poly_vec
        .iter()
        .zip(assignment.iter())
        .fold(Polynomial::<F>::zero(), |acc, (poly, &ref a)| {
            acc + poly.clone() * a.clone()
        })
}

fn get_h<F: Field>(qap: &QAP<F>, assignment: &Vec<F>) -> Polynomial<F> {
    let ell = accumulate_polynomials(&qap.ell_i_vec, assignment);
    let r = accumulate_polynomials(&qap.r_i_vec, assignment);
    let o = accumulate_polynomials(&qap.o_i_vec, assignment);
    (ell * r - o) / qap.t.clone()
}

pub fn setup(g1: &G1Point, g2: &G2Point, qap: &QAP<FqOrder>) -> (ProofKey, VerificationKey) {
    let s = FqOrder::random_element(&[]);
    let alpha = FqOrder::random_element(&[]);

    (
        ProofKey {
            g1_ell_i_vec: generate_challenge_vec(g1, &qap.ell_i_vec, &s),
            g1_r_i_vec: generate_challenge_vec(g1, &qap.r_i_vec, &s),
            g2_r_i_vec: generate_challenge_vec(g2, &qap.r_i_vec, &s),
            g1_o_i_vec: generate_challenge_vec(g1, &qap.o_i_vec, &s),
            g1_alpha_ell_i_vec: generate_alpha_challenge_vec(g1, &qap.ell_i_vec, &s, &alpha),
            g1_alpha_r_i_vec: generate_alpha_challenge_vec(g1, &qap.r_i_vec, &s, &alpha),
            g1_alpha_o_i_vec: generate_alpha_challenge_vec(g1, &qap.o_i_vec, &s, &alpha),
            g1_sj_vec: generate_s_powers(g1, &s, qap.m),
        },
        VerificationKey {
            g2_alpha: g2.mul_ref(alpha.get_value()),
            g2_t_s: g2.mul_ref(qap.t.eval(&s).sanitize().get_value()),
        },
    )
}

pub fn prove(assignment: &Vec<FqOrder>, proof_key: &ProofKey, qap: &QAP<FqOrder>) -> Proof {
    Proof {
        g1_ell: accumulate_curve_points(&proof_key.g1_ell_i_vec, assignment),
        g1_r: accumulate_curve_points(&proof_key.g1_r_i_vec, assignment),
        g2_r: accumulate_curve_points(&proof_key.g2_r_i_vec, assignment),
        g1_o: accumulate_curve_points(&proof_key.g1_o_i_vec, assignment),
        g1_ell_prime: accumulate_curve_points(&proof_key.g1_alpha_ell_i_vec, assignment),
        g1_r_prime: accumulate_curve_points(&proof_key.g1_alpha_r_i_vec, assignment),
        g1_o_prime: accumulate_curve_points(&proof_key.g1_alpha_o_i_vec, assignment),
        g1_h: get_h(qap, assignment).eval_with_powers_on_curve(&proof_key.g1_sj_vec),
    }
}

pub fn verify(
    g1: &G1Point,
    g2: &G2Point,
    proof: &Proof,
    verification_key: &VerificationKey,
) -> bool {
    let pairing1 = optimal_ate_pairing(&proof.g1_ell, &verification_key.g2_alpha);
    let pairing2 = optimal_ate_pairing(&proof.g1_ell_prime, &g2);
    if pairing1 != pairing2 {
        return false;
    }

    let pairing3 = optimal_ate_pairing(&proof.g1_r, &verification_key.g2_alpha);
    let pairing4 = optimal_ate_pairing(&proof.g1_r_prime, &g2);
    if pairing3 != pairing4 {
        return false;
    }

    let pairing5 = optimal_ate_pairing(&proof.g1_o, &verification_key.g2_alpha);
    let pairing6 = optimal_ate_pairing(&proof.g1_o_prime, &g2);
    if pairing5 != pairing6 {
        return false;
    }

    let pairing7 = optimal_ate_pairing(&proof.g1_ell, &proof.g2_r);
    let pairing8 = optimal_ate_pairing(&proof.g1_h, &verification_key.g2_t_s);
    let pairing9 = optimal_ate_pairing(&proof.g1_o, &g2);

    pairing7 == pairing8 * pairing9
}

pub fn interchange_attack(proof: &Proof) -> Proof {
    let mut new_proof = proof.clone();
    new_proof.g1_r = proof.g1_ell.clone();
    new_proof.g1_r_prime = proof.g1_ell_prime.clone();
    new_proof
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::modules::bn128::BN128;
    use crate::modules::field::{FiniteFieldElement, ModEIP197};
    use crate::modules::r1cs::R1CS;

    #[test]
    fn test_zksnark_naive_single_multiplication() {
        let left = vec![
            vec![
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::one(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
            ],
            vec![
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::one(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
            ],
            vec![
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::one(),
                FqOrder::zero(),
            ],
        ];

        let right = vec![
            vec![
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::one(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
            ],
            vec![
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::one(),
                FqOrder::zero(),
                FqOrder::zero(),
            ],
            vec![
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::one(),
            ],
        ];

        let out = vec![
            vec![
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::one(),
                FqOrder::zero(),
            ],
            vec![
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::one(),
            ],
            vec![
                FqOrder::zero(),
                FqOrder::one(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
                FqOrder::zero(),
            ],
        ];

        let v = vec![
            FqOrder::one(),
            FqOrder::from_value(210),
            FqOrder::from_value(2),
            FqOrder::from_value(3),
            FqOrder::from_value(5),
            FqOrder::from_value(7),
            FqOrder::from_value(6),
            FqOrder::from_value(35),
        ];
        let v_prime = vec![
            FqOrder::one(),
            FqOrder::from_value(211),
            FqOrder::from_value(2),
            FqOrder::from_value(3),
            FqOrder::from_value(5),
            FqOrder::from_value(7),
            FqOrder::from_value(6),
            FqOrder::from_value(35),
        ];

        let g1 = BN128::generator_g1();
        let g2 = BN128::generator_g2();

        let r1cs = R1CS::new(left, right, out);
        let qap = QAP::from_r1cs(&r1cs);

        let (proof_key, verification_key) = setup(&g1, &g2, &qap);

        let proof = prove(&v, &proof_key, &qap);
        assert!(verify(&g1, &g2, &proof, &verification_key));

        let proof_prime = prove(&v_prime, &proof_key, &qap);
        assert!(!verify(&g1, &g2, &proof_prime, &verification_key));

        let bogus_proof = interchange_attack(&proof);
        assert!(verify(&g1, &g2, &bogus_proof, &verification_key));
    }
}
