use crate::modules::algebra::curve::bn128::{optimal_ate_pairing, FqOrder, G1Point, G2Point};
use crate::modules::algebra::field::Field;
use crate::modules::algebra::ring::Ring;
use crate::modules::arithmetization::qap::QAP;
use crate::modules::zksnark::utils::{
    accumulate_curve_points, accumulate_polynomials, generate_alpha_challenge_vec,
    generate_challenge_vec, generate_s_powers, get_h,
};

#[derive(Debug, Clone)]
pub struct ProofKey2 {
    g1_ell_i_vec: Vec<G1Point>,
    g2_r_i_vec: Vec<G2Point>,
    g1_o_i_vec: Vec<G1Point>,
    g1_alpha_ell_i_vec: Vec<G1Point>,
    g2_alpha_r_i_vec: Vec<G2Point>,
    g1_alpha_o_i_vec: Vec<G1Point>,
    g1_sj_vec: Vec<G1Point>,
}

#[derive(Debug, Clone)]
pub struct VerificationKey2 {
    g2_alpha_ell: G2Point,
    g1_alpha_r: G1Point,
    g2_alpha_o: G2Point,
    g2_t_s: G2Point,
}

#[derive(Debug, Clone)]
pub struct Proof2 {
    g1_ell: G1Point,
    g2_r: G2Point,
    g1_o: G1Point,
    g1_ell_prime: G1Point,
    g2_r_prime: G2Point,
    g1_o_prime: G1Point,
    g1_h: G1Point,
}

pub fn setup(g1: &G1Point, g2: &G2Point, qap: &QAP<FqOrder>) -> (ProofKey2, VerificationKey2) {
    let s = FqOrder::random_element(&[]);
    let alpha_ell = FqOrder::random_element(&[]);
    let alpha_r = FqOrder::random_element(&[]);
    let alpha_o = FqOrder::random_element(&[]);

    (
        ProofKey2 {
            g1_ell_i_vec: generate_challenge_vec(g1, &qap.ell_i_vec, &s),
            g2_r_i_vec: generate_challenge_vec(g2, &qap.r_i_vec, &s),
            g1_o_i_vec: generate_challenge_vec(g1, &qap.o_i_vec, &s),
            g1_alpha_ell_i_vec: generate_alpha_challenge_vec(g1, &qap.ell_i_vec, &s, &alpha_ell),
            g2_alpha_r_i_vec: generate_alpha_challenge_vec(g2, &qap.r_i_vec, &s, &alpha_r),
            g1_alpha_o_i_vec: generate_alpha_challenge_vec(g1, &qap.o_i_vec, &s, &alpha_o),
            g1_sj_vec: generate_s_powers(g1, &s, qap.m),
        },
        VerificationKey2 {
            g2_alpha_ell: g2 * &alpha_ell,
            g1_alpha_r: g1 * &alpha_r,
            g2_alpha_o: g2 * &alpha_o,
            g2_t_s: g2 * qap.t.eval(&s).sanitize(),
        },
    )
}

pub fn prove(assignment: &Vec<FqOrder>, proof_key: &ProofKey2, qap: &QAP<FqOrder>) -> Proof2 {
    Proof2 {
        g1_ell: accumulate_curve_points(&proof_key.g1_ell_i_vec, assignment),
        g2_r: accumulate_curve_points(&proof_key.g2_r_i_vec, assignment),
        g1_o: accumulate_curve_points(&proof_key.g1_o_i_vec, assignment),
        g1_ell_prime: accumulate_curve_points(&proof_key.g1_alpha_ell_i_vec, assignment),
        g2_r_prime: accumulate_curve_points(&proof_key.g2_alpha_r_i_vec, assignment),
        g1_o_prime: accumulate_curve_points(&proof_key.g1_alpha_o_i_vec, assignment),
        g1_h: get_h(qap, assignment).eval_with_powers_on_curve(&proof_key.g1_sj_vec),
    }
}

pub fn verify(
    g1: &G1Point,
    g2: &G2Point,
    proof: &Proof2,
    verification_key: &VerificationKey2,
) -> bool {
    let pairing1 = optimal_ate_pairing(&proof.g1_ell, &verification_key.g2_alpha_ell);
    let pairing2 = optimal_ate_pairing(&proof.g1_ell_prime, &g2);
    if pairing1 != pairing2 {
        return false;
    }

    let pairing3 = optimal_ate_pairing(&verification_key.g1_alpha_r, &proof.g2_r);
    let pairing4 = optimal_ate_pairing(&g1, &proof.g2_r_prime);
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

    pairing7 == pairing8 * pairing9
}

pub fn inconsistent_variable_attack(
    assignment_ell: &Vec<FqOrder>,
    assignment_r: &Vec<FqOrder>,
    assignment_o: &Vec<FqOrder>,
    proof_key: &ProofKey2,
    qap: &QAP<FqOrder>,
) -> Proof2 {
    let ell = accumulate_polynomials(&qap.ell_i_vec, assignment_ell);
    let r = accumulate_polynomials(&qap.r_i_vec, assignment_r);
    let o = accumulate_polynomials(&qap.o_i_vec, assignment_o);
    let h = (ell * r - o) / qap.t.clone();

    Proof2 {
        g1_ell: accumulate_curve_points(&proof_key.g1_ell_i_vec, assignment_ell),
        g2_r: accumulate_curve_points(&proof_key.g2_r_i_vec, assignment_r),
        g1_o: accumulate_curve_points(&proof_key.g1_o_i_vec, assignment_o),
        g1_ell_prime: accumulate_curve_points(&proof_key.g1_alpha_ell_i_vec, assignment_ell),
        g2_r_prime: accumulate_curve_points(&proof_key.g2_alpha_r_i_vec, assignment_r),
        g1_o_prime: accumulate_curve_points(&proof_key.g1_alpha_o_i_vec, assignment_o),
        g1_h: h.eval_with_powers_on_curve(&proof_key.g1_sj_vec),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use num_traits::{One, Zero};

    use crate::modules::algebra::curve::bn128::BN128;
    use crate::modules::arithmetization::r1cs::R1CS;

    #[test]
    fn test_snark_ni_single_multiplication() {
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

        let wrong_proof = prove(&v_prime, &proof_key, &qap);
        assert!(!verify(&g1, &g2, &wrong_proof, &verification_key));

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

        let bogus_proof = inconsistent_variable_attack(&v_ell, &v_r, &v_o, &proof_key, &qap);
        assert!(verify(&g1, &g2, &bogus_proof, &verification_key));
    }
}
