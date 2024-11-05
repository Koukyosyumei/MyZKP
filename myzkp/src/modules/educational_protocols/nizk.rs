use crate::modules::curve::{weil_pairing, EllipticCurve, EllipticCurvePoint};
use crate::modules::field::Field;
use crate::modules::polynomial::Polynomial;
use num_bigint::BigInt;

pub struct ProofKey<F: Field> {
    alpha: Vec<EllipticCurvePoint<F>>,
    alpha_prime: Vec<EllipticCurvePoint<F>>,
}

pub struct VerificationKey<F: Field> {
    g_r: EllipticCurvePoint<F>,
    g_t_s: EllipticCurvePoint<F>,
}

pub struct Proof<F: Field> {
    u_prime: EllipticCurvePoint<F>,
    v_prime: EllipticCurvePoint<F>,
    w_prime: EllipticCurvePoint<F>,
}

pub struct TrustedSetup<F: Field> {
    proof_key: ProofKey<F>,
    verification_key: VerificationKey<F>,
}

impl<F: Field> TrustedSetup<F> {
    pub fn generate(
        curve: &EllipticCurve<F>,
        g: &EllipticCurvePoint<F>,
        t: &Polynomial<F>,
        n: usize,
        s: &F,
        r: &F,
    ) -> Self {
        let mut alpha = Vec::with_capacity(n);
        let mut alpha_prime = Vec::with_capacity(n);

        let mut s_power = F::one();
        for _ in 0..n {
            alpha.push(g.clone() * s_power.clone().get_value());
            alpha_prime.push(g.clone() * (s_power.clone() * r.clone()).get_value());
            s_power = s_power * s.clone();
        }

        let g_r = g.clone() * r.clone().get_value();
        let g_t_s = g.clone() * t.eval(s).get_value();

        TrustedSetup {
            proof_key: ProofKey { alpha, alpha_prime },
            verification_key: VerificationKey { g_r, g_t_s },
        }
    }
}

pub struct Prover<F: Field> {
    p: Polynomial<F>,
    h: Polynomial<F>,
    curve: EllipticCurve<F>,
    g: EllipticCurvePoint<F>,
}

impl<F: Field> Prover<F> {
    pub fn new(
        p: Polynomial<F>,
        t: Polynomial<F>,
        curve: EllipticCurve<F>,
        g: EllipticCurvePoint<F>,
    ) -> Self {
        let h = p.clone() / t;
        Prover { p, h, curve, g }
    }

    pub fn generate_proof(&self, proof_key: &ProofKey<F>, delta: &F) -> Proof<F> {
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

pub struct Verifier<F: Field> {
    curve: EllipticCurve<F>,
    g: EllipticCurvePoint<F>,
}

impl<F: Field> Verifier<F> {
    pub fn new(curve: EllipticCurve<F>, g: EllipticCurvePoint<F>) -> Self {
        Verifier { curve, g }
    }

    pub fn verify(&self, proof: &Proof<F>, vk: &VerificationKey<F>) -> bool {
        // Check e(u', g^r) = e(w', g)
        let pairing1 = weil_pairing(
            proof.u_prime.clone(),
            vk.g_r.clone(),
            BigInt::from(self.curve.order),
            Some(self.g.clone()),
        );
        let pairing2 = weil_pairing(
            proof.w_prime.clone(),
            self.g.clone(),
            BigInt::from(self.curve.order),
            Some(self.g.clone()),
        );
        let check1 = pairing1 == pairing2;

        // Check u' = v'^t
        let check2 = proof.u_prime == proof.v_prime * vk.g_t_s.clone().get_value();

        check1 && check2
    }
}

pub fn non_interactive_zkp_protocol<F: Field>(
    prover: &Prover<F>,
    verifier: &Verifier<F>,
    setup: &TrustedSetup<F>,
    delta: &F,
) -> bool {
    let proof = prover.generate_proof(&setup.proof_key, delta);
    verifier.verify(&proof, &setup.verification_key)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modules::field::{FiniteFieldElement, ModEIP197};
    use num_bigint::ToBigInt;

    #[test]
    fn test_non_interactive_zkp() {
        type F = FiniteFieldElement<ModEIP197>;

        // Set up the curve
        let a = F::from_value(30);
        let b = F::from_value(34);
        let curve = EllipticCurve { a, b };

        // Set up the generator point
        let g = EllipticCurvePoint::new(F::from_value(36), F::from_value(60), curve.clone());

        // Set up the polynomials
        let p = Polynomial::from_monomials(&[F::from_value(1), F::from_value(2), F::from_value(3)]);
        let t = Polynomial::from_monomials(&[F::from_value(1), F::from_value(2)]);

        // Generate trusted setup
        let s = F::from_value(123);
        let r = F::from_value(456);
        let setup = TrustedSetup::generate(&curve, &g, &t, 3, &s, &r);

        // Create prover and verifier
        let prover = Prover::new(p, t, curve.clone(), g.clone());
        let verifier = Verifier::new(curve, g);

        // Run the protocol
        let delta = F::from_value(789);
        let result = non_interactive_zkp_protocol(&prover, &verifier, &setup, &delta);

        assert!(result);
    }
}
