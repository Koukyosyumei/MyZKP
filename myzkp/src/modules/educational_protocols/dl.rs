use crate::modules::field::Field;
use crate::modules::polynomial::Polynomial;
use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_traits::Zero;
use std::str::FromStr;

pub struct Prover<F: Field> {
    pub p: Polynomial<F>,
    pub t: Polynomial<F>,
    pub h: Polynomial<F>,
}

pub struct Verifier<F: Field> {
    t: Polynomial<F>,
    s: F,
    g: F,
}

impl<F: Field> Prover<F> {
    pub fn new(p: Polynomial<F>, t: Polynomial<F>) -> Self {
        let h = p.clone() / t.clone();
        Prover { p, t, h }
    }

    pub fn compute_values(&self, alpha_powers: &[F]) -> (F, F) {
        let g_p = self.p.eval_with_powers(alpha_powers);
        let g_h = self.h.eval_with_powers(alpha_powers);
        (g_p, g_h)
    }
}

impl<F: Field> Verifier<F> {
    pub fn new(t: Polynomial<F>, generator: i128) -> Self {
        let mut rng = rand::thread_rng();
        let s = F::from_value(rng.gen_bigint_range(
            &BigInt::zero(),
            &BigInt::from_str("4835703278458516698824704").unwrap(), // 2^82
        ));
        let g = F::from_value(generator);
        Verifier { t, s, g }
    }

    pub fn generate_challenge(&self, max_degree: usize) -> Vec<F> {
        let mut alpha_powers = vec![];
        for i in 0..(max_degree + 1) {
            alpha_powers.push(
                self.g
                    .pow(self.s.clone().pow(i.to_bigint().unwrap()).get_value()),
            );
        }
        alpha_powers
    }

    pub fn verify(&self, u: &F, v: &F) -> bool {
        let t_s = self.t.eval(&self.s);
        u == &v.pow(t_s.get_value())
    }
}

// Simulating a malicious prover
pub struct MaliciousProver<F: Field> {
    t: Polynomial<F>,
}

impl<F: Field> MaliciousProver<F> {
    pub fn new(t: Polynomial<F>) -> Self {
        MaliciousProver { t }
    }

    pub fn compute_malicious_values(&self, alpha_powers: &[F]) -> (F, F) {
        let g_t = self.t.eval_with_powers(alpha_powers);
        let z = F::random_element(&[]);
        let g = &alpha_powers[0];
        let fake_v = g.pow(z.get_value());
        let fake_u = g_t.pow(z.get_value());
        (fake_u, fake_v)
    }
}

pub fn discrete_log_protocol<F: Field>(prover: &Prover<F>, verifier: &Verifier<F>) -> bool {
    // Step 1 & 2: Verifier generates a challenge
    let max_degree = prover.p.degree();
    let alpha_powers = verifier.generate_challenge(max_degree as usize);

    // Step 3: Prover computes and sends u = g^p and v = g^h
    let (u, v) = prover.compute_values(&alpha_powers);

    // Step 4: Verifier checks whether u = v^t
    verifier.verify(&u, &v)
}

pub fn malicious_discrete_log_protocol<F: Field>(
    prover: &MaliciousProver<F>,
    verifier: &Verifier<F>,
) -> bool {
    // Step 1 & 2: Verifier generates a challenge
    let max_degree = prover.t.degree() as usize;
    let alpha_powers = verifier.generate_challenge(max_degree as usize);

    // Step 3: Malicious Prover computes and sends fake u and v
    let (fake_u, fake_v) = prover.compute_malicious_values(&alpha_powers);

    // Step 4: Verifier checks whether u = v^t (which will pass for the fake values)
    verifier.verify(&fake_u, &fake_v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modules::field::{FiniteFieldElement, ModEIP197};
    use crate::modules::polynomial::Polynomial;

    #[test]
    fn test_dl_protocol() {
        const GENERATOR: i128 = 5; // A primitive root modulo MODULUS

        type F = FiniteFieldElement<ModEIP197>;

        // Create polynomials P(x) and T(x)
        let p =
            Polynomial::from_monomials(&[F::from_value(-1), F::from_value(-2), F::from_value(-3)]);
        let t = Polynomial::from_monomials(&[F::from_value(-1), F::from_value(-2)]);

        // Honest protocol
        let honest_prover = Prover::new(p.clone(), t.clone());
        let verifier = Verifier::new(t.clone(), GENERATOR);

        let honest_result = discrete_log_protocol(&honest_prover, &verifier);
        assert!(honest_result);

        // Malicious protocol
        let malicious_prover = MaliciousProver::new(t);
        let malicious_result = malicious_discrete_log_protocol(&malicious_prover, &verifier);
        assert!(malicious_result);
    }
}
