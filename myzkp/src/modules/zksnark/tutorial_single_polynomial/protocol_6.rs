use num_traits::One;

use crate::modules::algebra::curve::bn128::{optimal_ate_pairing, FqOrder, G1Point, G2Point};
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::ring::Ring;

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

pub fn setup(
    g1: &G1Point,
    g2: &G2Point,
    t: &Polynomial<FqOrder>,
    n: usize,
) -> (ProofKey, VerificationKey) {
    let s = FqOrder::random_element(&[]);
    let r = FqOrder::random_element(&[]);

    let mut alpha = Vec::with_capacity(n);
    let mut alpha_prime = Vec::with_capacity(n);

    let mut s_power = FqOrder::one();
    for _ in 0..1 + n {
        alpha.push(g1.mul_ref(s_power.clone().get_value()));
        alpha_prime.push(g1.mul_ref((s_power.clone() * r.clone()).get_value()));
        s_power = s_power * s.clone();
    }

    let g_r = g2.mul_ref(r.clone().get_value());
    let g_t_s = g2.mul_ref(t.eval(&s).get_value());

    (
        ProofKey { alpha, alpha_prime },
        VerificationKey { g_r, g_t_s },
    )
}

pub fn prove(p: &Polynomial<FqOrder>, t: &Polynomial<FqOrder>, proof_key: &ProofKey) -> Proof {
    let h = p.clone() / t.clone();
    let delta = FqOrder::random_element(&[]);

    let g_p = p.eval_with_powers_on_curve(&proof_key.alpha);
    let g_h = h.eval_with_powers_on_curve(&proof_key.alpha);
    let g_p_prime = p.eval_with_powers_on_curve(&proof_key.alpha_prime);

    Proof {
        u_prime: g_p * delta.get_value(),
        v_prime: g_h * delta.get_value(),
        w_prime: g_p_prime * delta.get_value(),
    }
}

pub fn verify(g2: &G2Point, proof: &Proof, vk: &VerificationKey) -> bool {
    // Check e(u', g^r) = e(w', g)
    let pairing1 = optimal_ate_pairing(&proof.u_prime, &vk.g_r);
    let pairing2 = optimal_ate_pairing(&proof.w_prime, g2);
    let check1 = pairing1 == pairing2;

    // Check e(u', g^t) = e(v', g)
    let pairing3 = optimal_ate_pairing(&proof.u_prime, g2);
    let pairing4 = optimal_ate_pairing(&proof.v_prime, &vk.g_t_s);
    let check2 = pairing3 == pairing4;

    check1 && check2
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::modules::algebra::curve::bn128::BN128;

    #[test]
    fn test_non_interactive_zkp() {
        let g1 = BN128::generator_g1();
        let g2 = BN128::generator_g2();

        // Set up the polynomials
        let p = Polynomial::from_monomials(&[
            FqOrder::from_value(-1),
            FqOrder::from_value(-2),
            FqOrder::from_value(-3),
        ]);
        let t = Polynomial::from_monomials(&[FqOrder::from_value(-1), FqOrder::from_value(-2)]);

        let (proof_key, verification_key) = setup(&g1, &g2, &t, 3);
        let proof = prove(&p, &t, &proof_key);
        assert!(verify(&g2, &proof, &verification_key));
    }
}
