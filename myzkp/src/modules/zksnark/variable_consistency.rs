use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_traits::One;
use num_traits::Zero;
use std::str::FromStr;

use crate::modules::bn128::{optimal_ate_pairing, Fq, Fq2, FqOrder, G1Point, G2Point};
use crate::modules::field::Field;
use crate::modules::polynomial::Polynomial;
use crate::modules::qap::QAP;
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
    g1_checksum_vec: Vec<G1Point>,
}

pub struct VerificationKey {
    g2_alpha_ell: G2Point,
    g2_alpha_r: G2Point,
    g2_alpha_o: G2Point,
    g2_beta_ell: G2Point,
    g2_beta_r: G2Point,
    g2_beta_o: G2Point,
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
    g1_z: G1Point,
}

pub fn setup(g1: &G1Point, g2: &G2Point, qap: &QAP<FqOrder>) -> (ProofKey, VerificationKey) {
    let mut rng = rand::thread_rng();
    let s = FqOrder::random_element(&[]);
    let alpha_ell = FqOrder::random_element(&[]);
    let alpha_r = FqOrder::random_element(&[]);
    let alpha_o = FqOrder::random_element(&[]);
    let beta_ell = FqOrder::random_element(&[]);
    let beta_r = FqOrder::random_element(&[]);
    let beta_o = FqOrder::random_element(&[]);

    let mut g1_ell_i_vec = Vec::with_capacity(qap.d);
    let mut g1_r_i_vec = Vec::with_capacity(qap.d);
    let mut g2_r_i_vec = Vec::with_capacity(qap.d);
    let mut g1_o_i_vec = Vec::with_capacity(qap.d);
    let mut g1_alpha_ell_i_vec = Vec::with_capacity(qap.d);
    let mut g1_alpha_r_i_vec = Vec::with_capacity(qap.d);
    let mut g1_alpha_o_i_vec = Vec::with_capacity(qap.d);
    let mut g1_sj_vec = Vec::with_capacity(qap.m);
    let mut g1_checksum_vec = Vec::with_capacity(qap.d);

    for i in 0..qap.d {
        let ell_i_s = qap.ell_i_vec[i].eval(&s).sanitize();
        let r_i_s = qap.r_i_vec[i].eval(&s).sanitize();
        let o_i_s = qap.o_i_vec[i].eval(&s).sanitize();

        g1_ell_i_vec.push(g1.mul_ref(ell_i_s.get_value()));
        g1_r_i_vec.push(g1.mul_ref(r_i_s.get_value()));
        g2_r_i_vec.push(g2.mul_ref(r_i_s.get_value()));
        g1_o_i_vec.push(g1.mul_ref(o_i_s.get_value()));
        g1_alpha_ell_i_vec.push(g1.mul_ref((alpha_ell.mul_ref(&ell_i_s)).get_value()));
        g1_alpha_r_i_vec.push(g1.mul_ref((alpha_r.mul_ref(&r_i_s)).get_value()));
        g1_alpha_o_i_vec.push(g1.mul_ref((alpha_o.mul_ref(&o_i_s)).get_value()));
        g1_checksum_vec.push(
            g1.mul_ref(
                ((beta_ell.mul_ref(&ell_i_s))
                    .add_ref(&beta_r.mul_ref(&r_i_s))
                    .add_ref(&beta_o.mul_ref(&o_i_s)))
                .get_value(),
            ),
        );
    }

    let mut s_power = FqOrder::one();
    for _ in 0..1 + qap.m {
        g1_sj_vec.push(g1.mul_ref(s_power.clone().get_value()));
        s_power = s_power * s.clone();
    }

    let g2_alpha_ell = g2.mul_ref(alpha_ell.get_value());
    let g2_alpha_r = g2.mul_ref(alpha_r.get_value());
    let g2_alpha_o = g2.mul_ref(alpha_o.get_value());
    let g2_beta_ell = g2.mul_ref(beta_ell.get_value());
    let g2_beta_r = g2.mul_ref(beta_r.get_value());
    let g2_beta_o = g2.mul_ref(beta_o.get_value());
    let g2_t_s = g2.mul_ref(qap.t.eval(&s).sanitize().get_value());

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
            g1_checksum_vec,
        },
        VerificationKey {
            g2_alpha_ell,
            g2_alpha_r,
            g2_alpha_o,
            g2_beta_ell,
            g2_beta_r,
            g2_beta_o,
            g2_t_s,
        },
    )
}

pub fn prove(
    g1: &G1Point,
    g2: &G2Point,
    assignment: &Vec<FqOrder>,
    proof_key: &ProofKey,
    qap: &QAP<FqOrder>,
) -> Proof {
    let mut g1_ell = G1Point::point_at_infinity();
    let mut g1_r = G1Point::point_at_infinity();
    let mut g2_r = G2Point::point_at_infinity();
    let mut g1_o = G1Point::point_at_infinity();
    let mut g1_ell_prime = G1Point::point_at_infinity();
    let mut g1_r_prime = G1Point::point_at_infinity();
    let mut g1_o_prime = G1Point::point_at_infinity();
    let mut g1_z = G1Point::point_at_infinity();

    for i in 0..qap.d {
        g1_ell = g1_ell + proof_key.g1_ell_i_vec[i].mul_ref(assignment[i].get_value());
        g1_r = g1_r + proof_key.g1_r_i_vec[i].mul_ref(assignment[i].get_value());
        g2_r = g2_r + proof_key.g2_r_i_vec[i].mul_ref(assignment[i].get_value());
        g1_o = g1_o + proof_key.g1_o_i_vec[i].mul_ref(assignment[i].get_value());
        g1_ell_prime =
            g1_ell_prime + proof_key.g1_alpha_ell_i_vec[i].mul_ref(assignment[i].get_value());
        g1_r_prime = g1_r_prime + proof_key.g1_alpha_r_i_vec[i].mul_ref(assignment[i].get_value());
        g1_o_prime = g1_o_prime + proof_key.g1_alpha_o_i_vec[i].mul_ref(assignment[i].get_value());
        g1_z = g1_z + proof_key.g1_checksum_vec[i].mul_ref(assignment[i].get_value());
    }

    let mut ell = Polynomial::<FqOrder>::zero();
    let mut r = Polynomial::<FqOrder>::zero();
    let mut o = Polynomial::<FqOrder>::zero();
    for i in 0..qap.d {
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
        g1_z,
    }
}

pub fn verify(
    g1: &G1Point,
    g2: &G2Point,
    proof: &Proof,
    verification_key: &VerificationKey,
) -> bool {
    let pairing1 = optimal_ate_pairing(&proof.g1_ell, &verification_key.g2_alpha_ell);
    let pairing2 = optimal_ate_pairing(&proof.g1_ell_prime, &g2);
    if pairing1 != pairing2 {
        return false;
    }

    let pairing3 = optimal_ate_pairing(&proof.g1_r, &verification_key.g2_alpha_r);
    let pairing4 = optimal_ate_pairing(&proof.g1_r_prime, &g2);
    if pairing3 != pairing4 {
        return false;
    }

    let pairing5 = optimal_ate_pairing(&proof.g1_o, &verification_key.g2_alpha_o);
    let pairing6 = optimal_ate_pairing(&proof.g1_o_prime, &g2);
    if pairing5 != pairing6 {
        return false;
    }

    let pairing7 = optimal_ate_pairing(&proof.g1_ell, &proof.g2_r);
    let pairing8 = optimal_ate_pairing(&proof.g1_h, &verification_key.g2_t_s);
    let pairing9 = optimal_ate_pairing(&proof.g1_o, &g2);

    if pairing7 != pairing8 * pairing9 {
        return false;
    }

    let pairing10 = optimal_ate_pairing(&proof.g1_ell, &verification_key.g2_beta_ell);
    let pairing11 = optimal_ate_pairing(&proof.g1_r, &verification_key.g2_beta_r);
    let pairing12 = optimal_ate_pairing(&proof.g1_o, &verification_key.g2_beta_o);
    let pairing13 = optimal_ate_pairing(&proof.g1_z, &g2);

    pairing10 * pairing11 * pairing12 == pairing13
}

pub fn modify_qap<F: Field>(qap: &QAP<F>) -> QAP<F> {
    QAP {
        m: qap.m,
        d: qap.d,
        ell_i_vec: qap.r_i_vec.clone(),
        r_i_vec: qap.ell_i_vec.clone(),
        o_i_vec: qap.o_i_vec.clone(),
        t: qap.t.clone(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::modules::bn128::BN128;
    use crate::modules::field::{FiniteFieldElement, ModEIP197};
    use crate::modules::r1cs::R1CS;

    #[test]
    fn test_zksnark_variable_consistency_single_multiplication() {
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

        let proof = prove(&g1, &g2, &v, &proof_key, &qap);
        assert!(verify(&g1, &g2, &proof, &verification_key));

        let proof_prime = prove(&&g1, &g2, &v_prime, &proof_key, &qap);
        assert!(!verify(&g1, &g2, &proof_prime, &verification_key));

        //let m_qap = modify_qap(&qap);
        //let bogus_proof = prove(&g1, &g2, &v, &proof_key, &verification_key, &m_qap);
        //assert!(!verify(&g1, &g2, &bogus_proof, &verification_key, &qap));
    }
}
