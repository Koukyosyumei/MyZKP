use num_bigint::ToBigInt;

use crate::modules::algebra::field::Field;
use crate::modules::algebra::polynomial::Polynomial;

pub struct Prover3<F: Field> {
    pub p: Polynomial<F>,
    pub t: Polynomial<F>,
    pub h: Polynomial<F>,
}

pub struct Verifier3<F: Field> {
    t: Polynomial<F>,
    s: F,
    g: F,
}

impl<F: Field> Prover3<F> {
    pub fn new(p: Polynomial<F>, t: Polynomial<F>) -> Self {
        let h = p.clone() / t.clone();
        Prover3 { p, t, h }
    }

    pub fn compute_values(&self, s_powers: &[F]) -> (F, F) {
        let g_p = self.p.eval_with_powers(s_powers);
        let g_h = self.h.eval_with_powers(s_powers);
        (g_p, g_h)
    }
}

impl<F: Field> Verifier3<F> {
    pub fn new(t: Polynomial<F>, generator: i128) -> Self {
        let s = F::random_element(&[]);
        let g = F::from_value(generator);
        Verifier3 { t, s, g }
    }

    pub fn generate_challenge(&self, max_degree: usize) -> Vec<F> {
        let mut s_powers = vec![];
        for i in 0..(max_degree + 1) {
            s_powers.push(
                self.g
                    .pow(self.s.clone().pow_m1(i.to_bigint().unwrap()).get_value()),
            );
        }
        s_powers
    }

    pub fn verify(&self, u: &F, v: &F) -> bool {
        let t_s = self.t.eval_m1(&self.s);
        u == &v.pow(t_s.get_value())
    }
}

// Simulating a malicious prover
pub struct MaliciousProver3<F: Field> {
    t: Polynomial<F>,
}

impl<F: Field> MaliciousProver3<F> {
    pub fn new(t: Polynomial<F>) -> Self {
        MaliciousProver3 { t }
    }

    pub fn compute_malicious_values(&self, s_powers: &[F]) -> (F, F) {
        let g_t = self.t.eval_with_powers(s_powers);
        let z = F::random_element(&[]);
        let g = &s_powers[0];
        let fake_v = g.pow(z.get_value());
        let fake_u = g_t.pow(z.get_value());
        (fake_u, fake_v)
    }
}

pub fn discrete_log_protocol<F: Field>(prover: &Prover3<F>, verifier: &Verifier3<F>) -> bool {
    // Step 1 & 2: Verifier3 generates a challenge
    let max_degree = prover.p.degree();
    let s_powers = verifier.generate_challenge(max_degree as usize);

    // Step 3: Prover3 computes and sends u = g^p and v = g^h
    let (u, v) = prover.compute_values(&s_powers);

    // Step 4: Verifier3 checks whether u = v^t
    verifier.verify(&u, &v)
}

pub fn malicious_discrete_log_protocol<F: Field>(
    prover: &MaliciousProver3<F>,
    verifier: &Verifier3<F>,
) -> bool {
    // Step 1 & 2: Verifier3 generates a challenge
    let max_degree = prover.t.degree() as usize;
    let s_powers = verifier.generate_challenge(max_degree as usize);

    // Step 3: Malicious Prover3 computes and sends fake u and v
    let (fake_u, fake_v) = prover.compute_malicious_values(&s_powers);

    // Step 4: Verifier3 checks whether u = v^t (which will pass for the fake values)
    verifier.verify(&fake_u, &fake_v)
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::modules::algebra::field::{FiniteFieldElement, ModEIP197};
    use crate::modules::algebra::polynomial::Polynomial;
    use crate::modules::algebra::ring::Ring;

    #[test]
    fn test_dl_protocol() {
        const GENERATOR: i128 = 5; // A primitive root modulo MODULUS

        type F = FiniteFieldElement<ModEIP197>;

        // Create polynomials P(x) and T(x)
        let p = Polynomial::from_monomials(&[F::from_value(1), F::from_value(2), F::from_value(3)]);
        let t = Polynomial::from_monomials(&[F::from_value(1), F::from_value(2)]);

        // Honest protocol
        let honest_prover = Prover3::new(p.clone(), t.clone());
        let verifier = Verifier3::new(t.clone(), GENERATOR);

        let honest_result = discrete_log_protocol(&honest_prover, &verifier);
        assert!(honest_result);

        // Malicious protocol
        let malicious_prover = MaliciousProver3::new(t);
        let malicious_result = malicious_discrete_log_protocol(&malicious_prover, &verifier);
        assert!(malicious_result);
    }
}
