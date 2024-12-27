use crate::modules::algebra::curve::bn128::{optimal_ate_pairing, FqOrder, G1Point, G2Point};
use crate::modules::algebra::field::Field;
use crate::modules::algebra::ring::Ring;
use crate::modules::arithmetization::qap::QAP;
use crate::modules::zksnark::utils::{
    accumulate_curve_points, accumulate_polynomials, generate_alpha_challenge_vec,
    generate_challenge_vec, generate_s_powers, get_h,
};

#[derive(Debug, Clone)]
pub struct ProofKey4 {
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

#[derive(Debug, Clone)]
pub struct VerificationKey4 {
    g2_alpha_ell: G2Point,
    g2_alpha_r: G2Point,
    g2_alpha_o: G2Point,
    g2_beta_ell_eta: G2Point,
    g2_beta_r_eta: G2Point,
    g2_beta_o_eta: G2Point,
    g2_t_s: G2Point,
    g2_eta: G2Point,
}

#[derive(Debug, Clone)]
pub struct Proof4 {
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

pub fn setup(g1: &G1Point, g2: &G2Point, qap: &QAP<FqOrder>) -> (ProofKey4, VerificationKey4) {
    let s = FqOrder::random_element(&[]);
    let alpha_ell = FqOrder::random_element(&[]);
    let alpha_r = FqOrder::random_element(&[]);
    let alpha_o = FqOrder::random_element(&[]);
    let beta_ell = FqOrder::random_element(&[]);
    let beta_r = FqOrder::random_element(&[]);
    let beta_o = FqOrder::random_element(&[]);
    let eta = FqOrder::random_element(&[]);

    let mut g1_checksum_vec = Vec::with_capacity(qap.d);

    for i in 0..qap.d {
        let ell_i_s = qap.ell_i_vec[i].eval(&s).sanitize();
        let r_i_s = qap.r_i_vec[i].eval(&s).sanitize();
        let o_i_s = qap.o_i_vec[i].eval(&s).sanitize();
        g1_checksum_vec.push(
            g1.mul_ref(
                ((beta_ell.mul_ref(&ell_i_s))
                    .add_ref(&beta_r.mul_ref(&r_i_s))
                    .add_ref(&beta_o.mul_ref(&o_i_s)))
                .get_value(),
            ),
        );
    }

    (
        ProofKey4 {
            g1_ell_i_vec: generate_challenge_vec(g1, &qap.ell_i_vec, &s),
            g1_r_i_vec: generate_challenge_vec(g1, &qap.r_i_vec, &s),
            g2_r_i_vec: generate_challenge_vec(g2, &qap.r_i_vec, &s),
            g1_o_i_vec: generate_challenge_vec(g1, &qap.o_i_vec, &s),
            g1_alpha_ell_i_vec: generate_alpha_challenge_vec(g1, &qap.ell_i_vec, &s, &alpha_ell),
            g1_alpha_r_i_vec: generate_alpha_challenge_vec(g1, &qap.r_i_vec, &s, &alpha_r),
            g1_alpha_o_i_vec: generate_alpha_challenge_vec(g1, &qap.o_i_vec, &s, &alpha_o),
            g1_sj_vec: generate_s_powers(g1, &s, qap.m),
            g1_checksum_vec: g1_checksum_vec,
        },
        VerificationKey4 {
            g2_alpha_ell: g2.mul_ref(alpha_ell.get_value()),
            g2_alpha_r: g2.mul_ref(alpha_r.get_value()),
            g2_alpha_o: g2.mul_ref(alpha_o.get_value()),
            g2_beta_ell_eta: g2.mul_ref(beta_ell.get_value()).mul_ref(eta.get_value()),
            g2_beta_r_eta: g2.mul_ref(beta_r.get_value()).mul_ref(eta.get_value()),
            g2_beta_o_eta: g2.mul_ref(beta_o.get_value()).mul_ref(eta.get_value()),
            g2_t_s: g2.mul_ref(qap.t.eval(&s).sanitize().get_value()),
            g2_eta: g2.mul_ref(eta.get_value()),
        },
    )
}

pub fn prove(assignment: &Vec<FqOrder>, proof_key: &ProofKey4, qap: &QAP<FqOrder>) -> Proof4 {
    Proof4 {
        g1_ell: accumulate_curve_points(&proof_key.g1_ell_i_vec, assignment),
        g1_r: accumulate_curve_points(&proof_key.g1_r_i_vec, assignment),
        g2_r: accumulate_curve_points(&proof_key.g2_r_i_vec, assignment),
        g1_o: accumulate_curve_points(&proof_key.g1_o_i_vec, assignment),
        g1_ell_prime: accumulate_curve_points(&proof_key.g1_alpha_ell_i_vec, assignment),
        g1_r_prime: accumulate_curve_points(&proof_key.g1_alpha_r_i_vec, assignment),
        g1_o_prime: accumulate_curve_points(&proof_key.g1_alpha_o_i_vec, assignment),
        g1_h: get_h(qap, assignment).eval_with_powers_on_curve(&proof_key.g1_sj_vec),
        g1_z: accumulate_curve_points(&proof_key.g1_checksum_vec, assignment),
    }
}

pub fn verify(g2: &G2Point, proof: &Proof4, verification_key: &VerificationKey4) -> bool {
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

    let pairing10 = optimal_ate_pairing(&proof.g1_ell, &verification_key.g2_beta_ell_eta);
    let pairing11 = optimal_ate_pairing(&proof.g1_r, &verification_key.g2_beta_r_eta);
    let pairing12 = optimal_ate_pairing(&proof.g1_o, &verification_key.g2_beta_o_eta);
    let pairing13 = optimal_ate_pairing(&proof.g1_z, &verification_key.g2_eta);

    pairing10 * pairing11 * pairing12 == pairing13
}

pub fn interchange_attack(proof: &Proof4) -> Proof4 {
    let mut new_proof = proof.clone();
    new_proof.g1_r = proof.g1_ell.clone();
    new_proof.g1_r_prime = proof.g1_ell_prime.clone();
    new_proof
}

pub fn inconsistent_variable_attack(
    assignment_ell: &Vec<FqOrder>,
    assignment_r: &Vec<FqOrder>,
    assignment_o: &Vec<FqOrder>,
    proof_key: &ProofKey4,
    qap: &QAP<FqOrder>,
) -> Proof4 {
    let ell = accumulate_polynomials(&qap.ell_i_vec, assignment_ell);
    let r = accumulate_polynomials(&qap.r_i_vec, assignment_r);
    let o = accumulate_polynomials(&qap.o_i_vec, assignment_o);
    let h = (ell * r - o) / qap.t.clone();

    Proof4 {
        g1_ell: accumulate_curve_points(&proof_key.g1_ell_i_vec, assignment_ell),
        g1_r: accumulate_curve_points(&proof_key.g1_r_i_vec, assignment_r),
        g2_r: accumulate_curve_points(&proof_key.g2_r_i_vec, assignment_r),
        g1_o: accumulate_curve_points(&proof_key.g1_o_i_vec, assignment_o),
        g1_ell_prime: accumulate_curve_points(&proof_key.g1_alpha_ell_i_vec, assignment_ell),
        g1_r_prime: accumulate_curve_points(&proof_key.g1_alpha_r_i_vec, assignment_r),
        g1_o_prime: accumulate_curve_points(&proof_key.g1_alpha_o_i_vec, assignment_o),
        g1_h: h.eval_with_powers_on_curve(&proof_key.g1_sj_vec),
        g1_z: accumulate_curve_points(&proof_key.g1_checksum_vec, assignment_ell),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use num_traits::{One, Zero};

    use crate::modules::algebra::curve::bn128::BN128;
    use crate::modules::arithmetization::r1cs::R1CS;

    #[test]
    fn test_snark_variable_consistency_single_multiplication() {
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
        assert!(verify(&g2, &proof, &verification_key));

        let proof_prime = prove(&v_prime, &proof_key, &qap);
        assert!(!verify(&g2, &proof_prime, &verification_key));

        let bogus_proof_1 = interchange_attack(&proof);
        assert!(!verify(&g2, &bogus_proof_1, &verification_key));

        let v_ell = vec![
            FqOrder::one(),
            FqOrder::from_value(210),
            FqOrder::from_value(2),
            FqOrder::from_value(3),
            FqOrder::from_value(5),
            FqOrder::from_value(7),
            FqOrder::from_value(6),
            FqOrder::from_value(35),
        ];

        let v_r = vec![
            FqOrder::one(),
            FqOrder::one(),
            FqOrder::one(),
            FqOrder::one(),
            FqOrder::one(),
            FqOrder::one(),
            FqOrder::one(),
            FqOrder::one(),
        ];

        let v_o = vec![
            FqOrder::one(),
            FqOrder::from_value(6),
            FqOrder::zero(),
            FqOrder::zero(),
            FqOrder::zero(),
            FqOrder::zero(),
            FqOrder::from_value(2),
            FqOrder::from_value(5),
        ];

        let bogus_proof_2 = inconsistent_variable_attack(&v_ell, &v_r, &v_o, &proof_key, &qap);
        assert!(!verify(&g2, &bogus_proof_2, &verification_key));
    }
}
