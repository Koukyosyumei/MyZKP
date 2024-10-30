use rand::Rng;

use crate::modules::field::{Field, FiniteFieldElement};
use crate::modules::polynomial::Polynomial;

struct Prover<F: Field> {
    p: Polynomial<F>,
    t: Polynomial<F>,
    h: Polynomial<F>,
}

struct Verifier<F: Field> {
    t: Polynomial<F>,
}

impl<F: Field> Prover<F> {
    fn new(p: Polynomial<F>, t: Polynomial<F>) -> Self {
        let h = p.clone() / t.clone();
        Prover { p, t, h }
    }

    fn compute_values(&self, s: &F) -> (F, F) {
        let h_s = self.h.eval(s);
        let p_s = self.p.eval(s);
        (h_s, p_s)
    }
}

impl<F: Field> Verifier<F> {
    fn new(t: Polynomial<F>) -> Self {
        Verifier { t }
    }

    fn generate_challenge(&self, modulus: i128) -> F {
        let mut rng = rand::thread_rng();
        F::from_value(rng.gen_range(0..modulus))
    }

    fn verify(&self, s: &F, h: &F, p: &F) -> bool {
        let t_s = self.t.eval(s);
        h.clone() * t_s == *p
    }
}

// Simulating a malicious prover
struct MaliciousProver<F: Field> {
    t: Polynomial<F>,
}

impl<F: Field> MaliciousProver<F> {
    fn new(t: Polynomial<F>) -> Self {
        MaliciousProver { t }
    }

    fn compute_malicious_values(&self, s: &F) -> (F, F) {
        let h_prime = F::random_element(&[]);
        let t_s = self.t.eval(s);
        let p_prime = h_prime.clone() * t_s;
        (h_prime, p_prime)
    }
}

fn schwartz_zippel_protocol<F: Field>(
    prover: &Prover<F>,
    verifier: &Verifier<F>,
    modulus: i128,
) -> bool {
    // Step 1: Verifier generates a random challenge
    let s = verifier.generate_challenge(modulus);

    // Step 2: Prover computes and sends h and p
    let (h, p) = prover.compute_values(&s);

    // Step 3: Verifier checks whether p = t * h
    verifier.verify(&s, &h, &p)
}

fn malicious_schwartz_zippel_protocol<F: Field>(
    prover: &MaliciousProver<F>,
    verifier: &Verifier<F>,
) -> bool {
    // Step 1: Verifier generates a random challenge
    let s = verifier.generate_challenge();

    // Step 2: Malicious Prover computes and sends h' and p'
    let (h_prime, p_prime) = prover.compute_malicious_values(&s);

    // Step 3: Verifier checks whether p' = t * h'
    verifier.verify(&s, &h_prime, &p_prime)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modules::field::FiniteFieldElement;
    use crate::modules::polynomial::Polynomial;

    #[test]
    fn test_sz_protocol() {
        const MODULUS: i128 = 3 * (1 << 30) + 1;
        type F = FiniteFieldElement<MODULUS>;

        // Create polynomials P(x) and T(x)
        let p = Polynomial::from_monomials(&[F::from_value(1), F::from_value(2), F::from_value(3)]);
        let t = Polynomial::from_monomials(&[F::from_value(1), F::from_value(2)]);

        // Honest protocol
        let honest_prover = Prover::new(p.clone(), t.clone());
        let verifier = Verifier::new(t.clone());

        let honest_result = schwartz_zippel_protocol(&honest_prover, &verifier);
        assert!(honest_result);

        // Malicious protocol
        let malicious_prover = MaliciousProver::new(t);
        let malicious_result = malicious_schwartz_zippel_protocol(&malicious_prover, &verifier);
        assert!(malicious_result);
    }
}
