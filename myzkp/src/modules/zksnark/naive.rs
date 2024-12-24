use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_traits::Zero;
use std::str::FromStr;

use crate::modules::bn128::{optimal_ate_pairing, Fq, G1Point, G2Point};
use crate::modules::field::Field;
use crate::modules::polynomial::Polynomial;
use crate::modules::ring::Ring;

pub struct ProofKey {
    g_ell_i_vec: Vec<G1Point>,
    g_r_i_vec: Vec<G1Point>,
    g_o_i_vec: Vec<G1Point>,
    g_alpha_ell_i_vec: Vec<G1Point>,
    g_alpha_r_i_vec: Vec<G1Point>,
    g_alpha_o_i_vec: Vec<G1Point>,
    g_sj_vec: Vec<G1Point>,
}

pub struct VerificationKey {
    g_alpha: G2Point,
    g_t_s: G2Point,
}

pub struct Proof {
    g_ell: G1Point,
    g_r: G1Point,
    g_o: G1Point,
    g_ell_prime: G1Point,
    g_r_prime: G1Point,
    g_o_prime: G1Point,
    g_h: G1Point,
}

pub struct QAP<F: Field> {
    pub m: usize,
    pub d: usize,
    pub ell_i_vec: Vec<Polynomial<F>>,
    pub r_i_vec: Vec<Polynomial<F>>,
    pub o_i_vec: Vec<Polynomial<F>>,
    pub t: Polynomial<F>,
}

pub fn setup<F: Field>(
    g1: &G1Point,
    g2: &G2Point,
    qap: &QAP<F>,
    n: usize,
) -> (ProofKey, VerificationKey) {
    let mut rng = rand::thread_rng();
    let s = Fq::from_value(rng.gen_bigint_range(&BigInt::zero(), &BigInt::from(std::u32::MAX)));
    let alpha = Fq::from_value(rng.gen_bigint_range(&BigInt::zero(), &BigInt::from(std::u32::MAX)));

    let mut g_ell_i_vec = Vec::with_capacity(qap.d);
    let mut g_r_i_vec = Vec::with_capacity(qap.d);
    let mut g_o_i_vec = Vec::with_capacity(qap.d);
    let mut g_alpha_ell_i_vec = Vec::with_capacity(qap.d);
    let mut g_alpha_r_i_vec = Vec::with_capacity(qap.d);
    let mut g_alpha_o_i_vec = Vec::with_capacity(qap.d);
    let mut g_sj_vec = Vec::with_capacity(n);

    for i in 0..1 + qap.d {
        g_ell_i_vec.push(g1.clone() * qap.ell_i_vec[i].eval(&s).get_value());
        g_r_i_vec.push(g1.clone() * qap.r_i_vec[i].eval(&s).get_value());
        g_o_i_vec.push(g1.clone() * qap.o_i_vec[i].eval(&s).get_value());
        g_alpha_ell_i_vec
            .push(g1.clone() * (alpha.clone() * qap.ell_i_vec[i].eval(&s).get_value()));
        g_alpha_r_i_vec.push(g1.clone() * (alpha.clone() * qap.r_i_vec[i].eval(&s).get_value()));
        g_alpha_o_i_vec.push(g1.clone() * (alpha.clone() * qap.o_i_vec[i].eval(&s).get_value()));
    }

    let mut s_power = Fq::one();
    for _ in 0..1 + n {
        g_sj_vec.push(g1.clone() * s_power.clone().get_value());
        s_power = s_power * s.clone();
    }

    let g_alpha = g2.clone() * alpha.clone().get_value();
    let g_t_s = g2.clone() * t.eval(&s).get_value();

    (
        ProofKey {
            g_ell_i_vec,
            g_r_i_vec,
            g_o_i_vec,
            g_alpha_ell_i_vec,
            g_alpha_r_i_vec,
            g_alpha_o_i_vec,
            g_sj_vec,
        },
        VerificationKey { g_alpha, g_t_s },
    )
}

pub fn generate_proof<F: Field>(
    &self,
    assignment: &Vec<F>,
    proof_key: &ProofKey,
    qap: &QAP<F>,
) -> Proof {
    let mut g_ell = proof_key.g_ell_i_vec[0] * assignment[0];
    let mut g_r = proof_key.g_r_i_vec[0] * assignment[0];
    let mut g_o = proof_key.g_o_i_vec[0] * assignment[0];
    let mut g_ell_prime = proof_key.g_alpha_ell_i_vec[0] * assignment[0];
    let mut g_r_prime = proof_key.g_alpha_r_i_vec[0] * assignment[0];
    let mut g_o_prime = proof_key.g_alpha_o_i_vec[0] * assignment[0];

    for i in 1..1 + qap.d {
        g_ell = g_ell + proof_key.g_ell_i_vec[i] * assignment[i];
        g_r = g_r + proof_key.g_r_i_vec[i] * assignment[i];
        g_o = g_o + proof_key.g_o_i_vec[i] * assignment[i];
        g_ell_prime = g_ell_prime + proof_key.g_alpha_ell_i_vec[i] * assignment[i];
        g_r_prime = g_r_prime + proof_key.g_alpha_r_i_vec[i] * assignment[i];
        g_o_prime = g_o_prime + proof_key.g_alpha_o_i_vec[i] * assignment[i];
    }

    let mut ell = Polynomial::<F>::zero();
    let mut r = Polynomial::<F>::zero();
    let mut o = Polynomial::<F>::zero();
    for i in 0..1 + qap.d {
        ell = ell + assignment[i] * qap.ell_i_vec[i];
        r = r + assignment[i] * qap.r_i_vec[i];
        o = o + assignment[i] * qap.o_i_vec[i];
    }
    let h = (rll * r - o) / qap.t;
}
