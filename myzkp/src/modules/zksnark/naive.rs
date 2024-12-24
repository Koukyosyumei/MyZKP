use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_traits::One;
use num_traits::Zero;
use std::str::FromStr;

use crate::modules::bn128::{optimal_ate_pairing, Fq, Fq2, G1Point, G2Point};
use crate::modules::field::Field;
use crate::modules::polynomial::Polynomial;
use crate::modules::ring::Ring;

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

pub struct VerificationKey {
    g2_alpha: G2Point,
    g2_t_s: G2Point,
}

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

pub struct QAP {
    pub m: usize,
    pub d: usize,
    pub ell_i_vec: Vec<Polynomial<Fq>>,
    pub r_i_vec: Vec<Polynomial<Fq>>,
    pub o_i_vec: Vec<Polynomial<Fq>>,
    pub t: Polynomial<Fq>,
}

pub fn setup(g1: &G1Point, g2: &G2Point, qap: &QAP, n: usize) -> (ProofKey, VerificationKey) {
    let mut rng = rand::thread_rng();
    let s = Fq::from_value(rng.gen_bigint_range(&BigInt::zero(), &BigInt::from(std::u32::MAX)));
    let alpha = Fq::from_value(rng.gen_bigint_range(&BigInt::zero(), &BigInt::from(std::u32::MAX)));

    let mut g1_ell_i_vec = Vec::with_capacity(qap.d);
    let mut g1_r_i_vec = Vec::with_capacity(qap.d);
    let mut g2_r_i_vec = Vec::with_capacity(qap.d);
    let mut g1_o_i_vec = Vec::with_capacity(qap.d);
    let mut g1_alpha_ell_i_vec = Vec::with_capacity(qap.d);
    let mut g1_alpha_r_i_vec = Vec::with_capacity(qap.d);
    let mut g1_alpha_o_i_vec = Vec::with_capacity(qap.d);
    let mut g1_sj_vec = Vec::with_capacity(n);

    for i in 0..1 + qap.d {
        g1_ell_i_vec.push(g1.clone() * qap.ell_i_vec[i].eval(&s).get_value());
        g1_r_i_vec.push(g1.clone() * qap.r_i_vec[i].eval(&s).get_value());
        g2_r_i_vec.push(g2.clone() * qap.r_i_vec[i].eval(&s).get_value());
        g1_o_i_vec.push(g1.clone() * qap.o_i_vec[i].eval(&s).get_value());
        g1_alpha_ell_i_vec
            .push(g1.clone() * (alpha.clone() * qap.ell_i_vec[i].eval(&s)).get_value());
        g1_alpha_r_i_vec.push(g1.clone() * (alpha.clone() * qap.r_i_vec[i].eval(&s)).get_value());
        g1_alpha_o_i_vec.push(g1.clone() * (alpha.clone() * qap.o_i_vec[i].eval(&s)).get_value());
    }

    let mut s_power = Fq::one();
    for _ in 0..1 + n {
        g1_sj_vec.push(g1.clone() * s_power.clone().get_value());
        s_power = s_power * s.clone();
    }

    let g2_alpha = g2.clone() * alpha.clone().get_value();
    let g2_t_s = g2.clone() * qap.t.eval(&s).get_value();

    (
        ProofKey {
            g1_ell_i_vec,
            g1_r_i_vec,
            g2_r_i_vec,
            g1_o_i_vec,
            g1_alpha_ell_i_vec,
            g1_alpha_r_i_vec,
            g1_alpha_o_i_vec,
            g1_sj_vec,
        },
        VerificationKey { g2_alpha, g2_t_s },
    )
}

pub fn prove(assignment: &Vec<Fq>, proof_key: &ProofKey, qap: &QAP) -> Proof {
    let mut g1_ell = proof_key.g1_ell_i_vec[0].clone() * assignment[0].get_value();
    let mut g1_r = proof_key.g1_r_i_vec[0].clone() * assignment[0].get_value();
    let mut g2_r = proof_key.g2_r_i_vec[0].clone() * assignment[0].get_value();
    let mut g1_o = proof_key.g1_o_i_vec[0].clone() * assignment[0].get_value();
    let mut g1_ell_prime = proof_key.g1_alpha_ell_i_vec[0].clone() * assignment[0].get_value();
    let mut g1_r_prime = proof_key.g1_alpha_r_i_vec[0].clone() * assignment[0].get_value();
    let mut g1_o_prime = proof_key.g1_alpha_o_i_vec[0].clone() * assignment[0].get_value();

    for i in 1..1 + qap.d {
        g1_ell = g1_ell + proof_key.g1_ell_i_vec[i].clone() * assignment[i].get_value();
        g1_r = g1_r + proof_key.g1_r_i_vec[i].clone() * assignment[i].get_value();
        g2_r = g2_r + proof_key.g2_r_i_vec[i].clone() * assignment[i].get_value();
        g1_o = g1_o + proof_key.g1_o_i_vec[i].clone() * assignment[i].get_value();
        g1_ell_prime =
            g1_ell_prime + proof_key.g1_alpha_ell_i_vec[i].clone() * assignment[i].get_value();
        g1_r_prime = g1_r_prime + proof_key.g1_alpha_r_i_vec[i].clone() * assignment[i].get_value();
        g1_o_prime = g1_o_prime + proof_key.g1_alpha_o_i_vec[i].clone() * assignment[i].get_value();
    }

    let mut ell = Polynomial::<Fq>::zero();
    let mut r = Polynomial::<Fq>::zero();
    let mut o = Polynomial::<Fq>::zero();
    for i in 0..1 + qap.d {
        ell = ell + qap.ell_i_vec[i].clone() * assignment[i].clone();
        r = r + qap.r_i_vec[i].clone() * assignment[i].clone();
        o = o + qap.o_i_vec[i].clone() * assignment[i].clone();
    }
    let h = (ell * r - o) / qap.t.clone();
    let g1_h = h.eval_with_powers_on_curve(&proof_key.g1_sj_vec);

    Proof {
        g1_ell,
        g1_r,
        g2_r,
        g1_o,
        g1_ell_prime,
        g1_r_prime,
        g1_o_prime,
        g1_h,
    }
}

pub fn verify(
    g1: &G1Point,
    g2: &G2Point,
    proof: &Proof,
    verification_key: &VerificationKey,
    qap: &QAP,
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
