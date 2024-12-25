use num_bigint::{BigInt, RandBigInt};
use num_traits::One;
use num_traits::Zero;
use std::str::FromStr;

use crate::modules::bn128::{optimal_ate_pairing, Fq, G1Point, G2Point};
use crate::modules::polynomial::Polynomial;
use crate::modules::ring::Ring;

pub struct ProofKey {
    alpha: Vec<G1Point>,
    alpha_prime: Vec<G1Point>,
}

pub struct VerificationKey {
    g_r: G2Point,
    g_t_s: G2Point,
}

pub struct Proof {
    u_prime: G1Point,
    v_prime: G1Point,
    w_prime: G1Point,
}

pub struct TrustedSetup {
    proof_key: ProofKey,
    verification_key: VerificationKey,
}

impl TrustedSetup {
    pub fn generate(g1: &G1Point, g2: &G2Point, t: &Polynomial<Fq>, n: usize) -> Self {
        let mut rng = rand::thread_rng();
        let s = Fq::from_value(rng.gen_bigint_range(&BigInt::zero(), &BigInt::from(std::u32::MAX)));
        let r = Fq::from_value(rng.gen_bigint_range(&BigInt::zero(), &BigInt::from(std::u32::MAX)));

        let mut alpha = Vec::with_capacity(n);
        let mut alpha_prime = Vec::with_capacity(n);

        let mut s_power = Fq::one();
        for _ in 0..1 + n {
            alpha.push(g1.mul_ref(s_power.clone().get_value()));
            alpha_prime.push(g1.mul_ref((s_power.clone() * r.clone()).get_value()));
            s_power = s_power * s.clone();
        }

        let g_r = g2.mul_ref(r.clone().get_value());
        let g_t_s = g2.mul_ref(t.eval(&s).get_value());

        TrustedSetup {
            proof_key: ProofKey { alpha, alpha_prime },
            verification_key: VerificationKey { g_r, g_t_s },
        }
    }
}

pub struct Prover {
    pub p: Polynomial<Fq>,
    pub h: Polynomial<Fq>,
}

impl Prover {
    pub fn new(p: Polynomial<Fq>, t: Polynomial<Fq>) -> Self {
        let h = p.clone() / t;
        Prover { p, h }
    }

    pub fn generate_proof(&self, proof_key: &ProofKey) -> Proof {
        let mut rng = rand::thread_rng();
        let delta =
            Fq::from_value(rng.gen_bigint_range(&BigInt::zero(), &BigInt::from(std::u32::MAX)));

        let g_p = self.p.eval_with_powers_on_curve(&proof_key.alpha);
        let g_h = self.h.eval_with_powers_on_curve(&proof_key.alpha);
        let g_p_prime = self.p.eval_with_powers_on_curve(&proof_key.alpha_prime);

        Proof {
            u_prime: g_p * delta.get_value(),
            v_prime: g_h * delta.get_value(),
            w_prime: g_p_prime * delta.get_value(),
        }
    }
}

pub struct Verifier {
    pub g1: G1Point,
    pub g2: G2Point,
}

impl Verifier {
    pub fn new(g1: G1Point, g2: G2Point) -> Self {
        Verifier { g1, g2 }
    }

    pub fn verify(&self, proof: &Proof, vk: &VerificationKey) -> bool {
        // Check e(u', g^r) = e(w', g)
        let pairing1 = optimal_ate_pairing(&proof.u_prime, &vk.g_r);
        let pairing2 = optimal_ate_pairing(&proof.w_prime, &self.g2);
        let check1 = pairing1 == pairing2;

        // Check e(u', g^t) = e(v', g)
        let pairing3 = optimal_ate_pairing(&proof.u_prime, &self.g2);
        let pairing4 = optimal_ate_pairing(&proof.v_prime, &vk.g_t_s);
        let check2 = pairing3 == pairing4;

        check1 && check2
    }
}

pub fn non_interactive_zkp_protocol(
    prover: &Prover,
    verifier: &Verifier,
    setup: &TrustedSetup,
) -> bool {
    let proof = prover.generate_proof(&setup.proof_key);
    verifier.verify(&proof, &setup.verification_key)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modules::bn128::BN128;

    #[test]
    fn test_non_interactive_zkp() {
        let g1 = BN128::generator_g1();
        let g2 = BN128::generator_g2();

        // Set up the polynomials
        let p = Polynomial::from_monomials(&[
            Fq::from_value(-1),
            Fq::from_value(-2),
            Fq::from_value(-3),
        ]);
        let t = Polynomial::from_monomials(&[Fq::from_value(-1), Fq::from_value(-2)]);

        let setup = TrustedSetup::generate(&g1, &g2, &t, 3);

        // Create prover and verifier
        let prover = Prover::new(p, t.clone());
        let verifier = Verifier::new(g1.clone(), g2.clone());

        // Run the protocol
        let result = non_interactive_zkp_protocol(&prover, &verifier, &setup);
        assert!(result);
    }
}
