use crate::modules::curve::{weil_pairing, EllipticCurve, EllipticCurvePoint};
use crate::modules::field::Field;
use crate::modules::polynomial::Polynomial;
use num_bigint::BigInt;

pub struct ProofKey<F: Field, E: EllipticCurve> {
    alpha: Vec<EllipticCurvePoint<F, E>>,
    alpha_prime: Vec<EllipticCurvePoint<F, E>>,
}

pub struct VerificationKey<F: Field, E: EllipticCurve> {
    g_r: EllipticCurvePoint<F, E>,
    g_t_s: EllipticCurvePoint<F, E>,
}

pub struct Proof<F: Field, E: EllipticCurve> {
    u_prime: EllipticCurvePoint<F, E>,
    v_prime: EllipticCurvePoint<F, E>,
    w_prime: EllipticCurvePoint<F, E>,
}

pub struct TrustedSetup<F: Field, E: EllipticCurve> {
    proof_key: ProofKey<F, E>,
    verification_key: VerificationKey<F, E>,
}

impl<F: Field, E: EllipticCurve> TrustedSetup<F, E> {
    pub fn generate(
        g: &EllipticCurvePoint<F, E>,
        t: &Polynomial<F>,
        n: usize,
        s: &F,
        r: &F,
    ) -> Self {
        let mut alpha = Vec::with_capacity(n);
        let mut alpha_prime = Vec::with_capacity(n);

        let mut s_power = F::one();
        for _ in 0..1 + n {
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

pub struct Prover<F: Field, E: EllipticCurve> {
    pub p: Polynomial<F>,
    pub h: Polynomial<F>,
    pub g: EllipticCurvePoint<F, E>,
}

impl<F: Field, E: EllipticCurve> Prover<F, E> {
    pub fn new(p: Polynomial<F>, t: Polynomial<F>, g: EllipticCurvePoint<F, E>) -> Self {
        let h = p.clone() / t;
        Prover { p, h, g }
    }

    pub fn generate_proof(&self, proof_key: &ProofKey<F, E>, delta: &F) -> Proof<F, E> {
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

pub struct Verifier<F: Field, E: EllipticCurve> {
    g: EllipticCurvePoint<F, E>,
}

impl<F: Field, E: EllipticCurve> Verifier<F, E> {
    pub fn new(g: EllipticCurvePoint<F, E>) -> Self {
        Verifier { g }
    }

    pub fn verify(
        &self,
        proof: &Proof<F, E>,
        vk: &VerificationKey<F, E>,
        curve_order: BigInt,
        s: EllipticCurvePoint<F, E>,
    ) -> bool {
        // Check e(u', g^r) = e(w', g)
        let pairing1 = weil_pairing(
            proof.u_prime.clone(),
            vk.g_r.clone(),
            BigInt::from(curve_order.clone()),
            Some(self.g.clone()),
        );
        let pairing2 = weil_pairing(
            proof.w_prime.clone(),
            self.g.clone(),
            BigInt::from(curve_order.clone()),
            Some(self.g.clone()),
        );
        let check1 = pairing1 == pairing2;

        println!("^^^^ {} {} {}", pairing1, pairing2, check1);

        println!(
            "{}",
            proof.u_prime.clone() == self.g.clone() * BigInt::from(24)
        );

        println!(
            "{}",
            proof.v_prime.clone() == self.g.clone() * BigInt::from(4)
        );

        println!("{}", vk.g_t_s.clone() == self.g.clone() * BigInt::from(6));

        // Check e(u', g^t) = e(v', g)
        let pairing3 = weil_pairing(
            proof.u_prime.clone(),
            self.g.clone(),
            BigInt::from(curve_order.clone()),
            Some(s.clone()),
        );
        let pairing4 = weil_pairing(
            proof.v_prime.clone(),
            vk.g_t_s.clone(),
            BigInt::from(curve_order.clone()),
            Some(s.clone()),
        );
        let check2 = pairing3 == pairing4;

        println!("^^^^ {} {} {}", pairing3, pairing4, check2);

        check1 && check2
    }
}

pub fn non_interactive_zkp_protocol<F: Field, E: EllipticCurve>(
    prover: &Prover<F, E>,
    verifier: &Verifier<F, E>,
    setup: &TrustedSetup<F, E>,
    delta: &F,
    curve_order: BigInt,
    ms: EllipticCurvePoint<F, E>,
) -> bool {
    let proof = prover.generate_proof(&setup.proof_key, delta);
    verifier.verify(&proof, &setup.verification_key, curve_order, ms)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    use crate::modules::field::ModulusValue;
    use crate::modules::ring::Ring;
    use crate::{
        define_myzkp_curve_type, define_myzkp_modulus_type,
        modules::field::{FiniteFieldElement, ModEIP197},
    };

    #[test]
    fn test_non_interactive_zkp() {
        define_myzkp_modulus_type!(Mod631, "631");
        define_myzkp_curve_type!(BLS12to381, "1", "4");
        type F = FiniteFieldElement<Mod631>;
        type E = CurveA30B34;

        // Set up the generator point
        let g = EllipticCurvePoint::<F, E>::new(F::from_value(36), F::from_value(60));

        // Set up the polynomials
        let p =
            Polynomial::from_monomials(&[F::from_value(-1), F::from_value(-2), F::from_value(-3)]);
        let t = Polynomial::from_monomials(&[F::from_value(-1), F::from_value(-2)]);

        // Generate trusted setup
        let s = F::from_value(2);
        let r = F::from_value(1); //F::from_value(456);
        let setup = TrustedSetup::generate(&g, &t, 3, &s, &r);

        // Create prover and verifier
        let prover = Prover::new(p, t.clone(), g.clone());
        let verifier = Verifier::new(g.clone());

        let ms = EllipticCurvePoint::<F, E>::new(F::from_value(0_i64), F::from_value(36_i64));

        // Run the protocol
        let delta = F::from_value(23);
        let order = BigInt::from_str("5").unwrap();
        let result = non_interactive_zkp_protocol(&prover, &verifier, &setup, &delta, order, ms);

        let proof = prover.generate_proof(&setup.proof_key, &delta);
        /*
        println!(
            "v_prime={}, {}",
            proof.v_prime.clone(),
            proof.v_prime == g.clone() * prover.h.eval(&s).get_value()
        );

        println!(
            "p(s)={} h(s)={} t(s)={}",
            prover.p.eval(&s).get_value(),
            prover.h.eval(&s).get_value(),
            t.eval(&s).get_value()
        );
        */
        println!(
            "@@  {}",
            //proof.u_prime.clone() * r.get_value(),
            //proof.w_prime.clone(),
            proof.u_prime.clone() * r.get_value() == proof.w_prime
        );
        println!(
            "@@  {}",
            //proof.u_prime.clone(),
            //proof.v_prime.clone() * t.eval(&s).get_value(),
            proof.u_prime.clone() == proof.v_prime * t.eval(&s).get_value()
        );

        assert!(result);
    }
}
