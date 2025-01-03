use num_traits::One;

use crate::modules::algebra::curve::bn128::{optimal_ate_pairing, FqOrder, G1Point, G2Point};
use crate::modules::algebra::field::Field;
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::ring::Ring;
use crate::modules::arithmetization::qap::QAP;
use crate::modules::zksnark::utils::{
    accumulate_curve_points, accumulate_polynomials, generate_alpha_challenge_vec,
    generate_challenge_vec, generate_s_powers,
};

pub struct PinocchioProofKey {
    g1_ell_i_vec: Vec<G1Point>,
    g2_r_i_vec: Vec<G2Point>,
    g1_o_i_vec: Vec<G1Point>,
    g1_alpha_ell_i_vec: Vec<G1Point>,
    g2_alpha_r_i_vec: Vec<G2Point>,
    g1_alpha_o_i_vec: Vec<G1Point>,
    g1_sj_vec: Vec<G1Point>,
    g1_checksum_vec: Vec<G1Point>,
    g1_ell_ts: G1Point,
    g2_r_ts: G2Point,
    g1_o_ts: G1Point,
    g1_ell_alpha_ts: G1Point,
    g2_r_alpha_ts: G2Point,
    g1_o_alpha_ts: G1Point,
    g1_ell_beta_ts: G1Point,
    g1_r_beta_ts: G1Point,
    g1_o_beta_ts: G1Point,
}

pub struct PinocchioVerificationKey {
    g2_alpha_ell: G2Point,
    g1_alpha_r: G1Point,
    g2_alpha_o: G2Point,
    g1_beta_eta: G1Point,
    g2_beta_eta: G2Point,
    g2_t_s: G2Point,
    g2_eta: G2Point,
}

pub struct PinocchioProof {
    g1_ell: G1Point,
    g2_r: G2Point,
    g1_o: G1Point,
    g1_ell_prime: G1Point,
    g2_r_prime: G2Point,
    g1_o_prime: G1Point,
    g1_h: G1Point,
    g1_z: G1Point,
}

/// Generates the proving and verification keys for the Pinocchio zk-SNARK scheme.
///
/// # Parameters:
/// - `g1`: A reference to a G1 curve point.
/// - `g2`: A reference to a G2 curve point.
/// - `qap`: The quadratic arithmetic program (QAP) representing the circuit.
///
/// # Returns:
/// A tuple containing:
/// - `PinocchioProofKey`: The proving key.
/// - `PinocchioVerificationKey`: The verification key.
pub fn setup(
    g1: &G1Point,
    g2: &G2Point,
    qap: &QAP<FqOrder>,
) -> (PinocchioProofKey, PinocchioVerificationKey) {
    let s = FqOrder::random_element(&[]);
    let alpha_ell = FqOrder::random_element(&[]);
    let alpha_r = FqOrder::random_element(&[]);
    let alpha_o = FqOrder::random_element(&[]);
    let beta = FqOrder::random_element(&[]);
    let eta = FqOrder::random_element(&[]);
    let rho_ell = FqOrder::random_element(&[]);
    let rho_r = FqOrder::random_element(&[]);
    let rho_o = &rho_ell * &rho_r;

    let g1_ell = g1 * &rho_ell;
    let g1_r = g1 * &rho_r;
    let g2_r = g2 * &rho_r;
    let g1_o = g1 * &rho_o;
    let g2_o = g2 * &rho_o;

    let mut g1_checksum_vec = Vec::with_capacity(qap.d);

    for i in 0..qap.d {
        let ell_i_s = qap.ell_i_vec[i].eval(&s).sanitize();
        let r_i_s = qap.r_i_vec[i].eval(&s).sanitize();
        let o_i_s = qap.o_i_vec[i].eval(&s).sanitize();
        g1_checksum_vec.push(
            &g1_ell * (&beta * &ell_i_s) + &g1_r * (&beta * &r_i_s) + &g1_o * (&beta * &o_i_s),
        );
    }

    let t_s = qap.t.eval(&s).sanitize();
    let beta_eta = &beta * &eta;
    let beta_t_s = beta * &t_s;

    (
        PinocchioProofKey {
            g1_ell_i_vec: generate_challenge_vec(&g1_ell, &qap.ell_i_vec, &s),
            g2_r_i_vec: generate_challenge_vec(&g2_r, &qap.r_i_vec, &s),
            g1_o_i_vec: generate_challenge_vec(&g1_o, &qap.o_i_vec, &s),
            g1_alpha_ell_i_vec: generate_alpha_challenge_vec(
                &g1_ell,
                &qap.ell_i_vec,
                &s,
                &alpha_ell,
            ),
            g2_alpha_r_i_vec: generate_alpha_challenge_vec(&g2_r, &qap.r_i_vec, &s, &alpha_r),
            g1_alpha_o_i_vec: generate_alpha_challenge_vec(&g1_o, &qap.o_i_vec, &s, &alpha_o),
            g1_sj_vec: generate_s_powers(&g1, &s, qap.m),
            g1_checksum_vec: g1_checksum_vec,
            g1_ell_ts: &g1_ell * &t_s,
            g2_r_ts: &g2_r * &t_s,
            g1_o_ts: &g1_o * &t_s,
            g1_ell_alpha_ts: &g1_ell * (&t_s * &alpha_ell),
            g2_r_alpha_ts: &g2_r * (&t_s * &alpha_r),
            g1_o_alpha_ts: &g1_o * (&t_s * &alpha_o),
            g1_ell_beta_ts: g1_ell * beta_t_s.get_value(),
            g1_r_beta_ts: g1_r * beta_t_s.get_value(),
            g1_o_beta_ts: g1_o * beta_t_s.get_value(),
        },
        PinocchioVerificationKey {
            g2_alpha_ell: g2 * &alpha_ell,
            g1_alpha_r: g1 * &alpha_r,
            g2_alpha_o: g2 * &alpha_o,
            g1_beta_eta: g1 * &beta_eta,
            g2_beta_eta: g2 * &beta_eta,
            g2_t_s: g2_o * t_s.get_value(),
            g2_eta: g2 * &eta,
        },
    )
}

pub fn get_shifted_h(
    qap: &QAP<FqOrder>,
    assignment: &Vec<FqOrder>,
    delta_ell: &FqOrder,
    delta_r: &FqOrder,
    delta_o: &FqOrder,
) -> Polynomial<FqOrder> {
    let ell = accumulate_polynomials(&qap.ell_i_vec, assignment);
    let r = accumulate_polynomials(&qap.r_i_vec, assignment);
    let o = accumulate_polynomials(&qap.o_i_vec, assignment);
    (&ell * &r - o) / qap.t.clone()
        + ell * delta_r
        + r * delta_ell
        + qap.t.clone() * (delta_ell * delta_r)
        - Polynomial::<FqOrder>::one() * delta_o
}

/// Generates a zk-SNARK proof using the Pinocchio scheme.
///
/// # Parameters:
/// - `assignment`: A vector of field elements representing the circuit assignment.
/// - `proof_key`: A reference to the proving key.
/// - `qap`: The QAP representation of the circuit.
///
/// # Returns:
/// A `PinocchioProof` object representing the proof.
pub fn prove(
    assignment: &Vec<FqOrder>,
    proof_key: &PinocchioProofKey,
    qap: &QAP<FqOrder>,
) -> PinocchioProof {
    let delta_ell = FqOrder::random_element(&[]);
    let delta_r = FqOrder::random_element(&[]);
    let delta_o = FqOrder::random_element(&[]);

    PinocchioProof {
        g1_ell: &proof_key.g1_ell_ts * &delta_ell
            + accumulate_curve_points(&proof_key.g1_ell_i_vec, assignment),
        g2_r: &proof_key.g2_r_ts * &delta_r
            + accumulate_curve_points(&proof_key.g2_r_i_vec, assignment),
        g1_o: &proof_key.g1_o_ts * &delta_o
            + accumulate_curve_points(&proof_key.g1_o_i_vec, assignment),
        g1_ell_prime: &proof_key.g1_ell_alpha_ts * &delta_ell
            + accumulate_curve_points(&proof_key.g1_alpha_ell_i_vec, assignment),
        g2_r_prime: &proof_key.g2_r_alpha_ts * &delta_r
            + accumulate_curve_points(&proof_key.g2_alpha_r_i_vec, assignment),
        g1_o_prime: &proof_key.g1_o_alpha_ts * &delta_o
            + accumulate_curve_points(&proof_key.g1_alpha_o_i_vec, assignment),
        g1_h: get_shifted_h(qap, assignment, &delta_ell, &delta_r, &delta_o)
            .eval_with_powers_on_curve(&proof_key.g1_sj_vec),
        g1_z: &proof_key.g1_ell_beta_ts * &delta_ell
            + &proof_key.g1_r_beta_ts * &delta_r
            + &proof_key.g1_o_beta_ts * &delta_o
            + accumulate_curve_points(&proof_key.g1_checksum_vec, assignment),
    }
}

/// Verifies a zk-SNARK proof using the Pinocchio scheme.
///
/// # Parameters:
/// - `g1`: A reference to a G1 curve point.
/// - `g2`: A reference to a G2 curve point.
/// - `proof`: A reference to the proof to be verified.
/// - `verification_key`: A reference to the verification key.
///
/// # Returns:
/// `true` if the proof is valid; otherwise, `false`.
pub fn verify(
    g1: &G1Point,
    g2: &G2Point,
    proof: &PinocchioProof,
    verification_key: &PinocchioVerificationKey,
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

    if pairing7 != pairing8 * pairing9 {
        return false;
    }

    let pairing10 = optimal_ate_pairing(
        &proof.g1_ell.add_ref(&proof.g1_o),
        &verification_key.g2_beta_eta,
    );
    let pairing11 = optimal_ate_pairing(&verification_key.g1_beta_eta, &proof.g2_r);
    let pairing12 = optimal_ate_pairing(&proof.g1_z, &verification_key.g2_eta);

    pairing10 * pairing11 == pairing12
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
        assert!(verify(&g1, &g2, &proof, &verification_key));

        let wrong_proof = prove(&v_prime, &proof_key, &qap);
        assert!(!verify(&g1, &g2, &wrong_proof, &verification_key));
    }
}
