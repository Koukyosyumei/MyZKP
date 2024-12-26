use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_traits::Zero;
use std::str::FromStr;

use crate::modules::field::Field;
use crate::modules::polynomial::Polynomial;

pub struct Prover<F: Field> {
    pub p: Polynomial<F>,
    pub t: Polynomial<F>,
    pub h: Polynomial<F>,
}

pub struct Verifier<F: Field> {
    t: Polynomial<F>,
    s: F,
    r: F,
    g: F,
}

impl<F: Field> Prover<F> {
    pub fn new(p: Polynomial<F>, t: Polynomial<F>) -> Self {
        let h = p.clone() / t.clone();
        Prover { p, t, h }
    }

    pub fn compute_values(&self, alpha_powers: &[F], alpha_prime_powers: &[F]) -> (F, F, F) {
        let g_p = self.p.eval_with_powers(alpha_powers);
        let g_h = self.h.eval_with_powers(alpha_powers);
        let g_p_prime = self.p.eval_with_powers(alpha_prime_powers);
        (g_p, g_h, g_p_prime)
    }
}

impl<F: Field> Verifier<F> {
    pub fn new(t: Polynomial<F>, generator: i128) -> Self {
        let s = F::random_element(&[]);
        let r = F::random_element(&[]);
        let g = F::from_value(generator);
        Verifier { t, s, r, g }
    }

    pub fn generate_challenge(&self, max_degree: usize) -> (Vec<F>, Vec<F>) {
        let mut alpha_powers = vec![];
        let mut alpha_prime_powers = vec![];

        for i in 0..(max_degree + 1) {
            alpha_powers.push(
                self.g
                    .pow(self.s.clone().pow_m1(i.to_bigint().unwrap()).get_value()),
            );
            alpha_prime_powers.push(alpha_powers.last().unwrap().pow(self.r.get_value()));
        }

        (alpha_powers, alpha_prime_powers)
    }

    pub fn verify(&self, u: &F, v: &F, w: &F) -> bool {
        let t_s = self.t.eval_m1(&self.s);
        let u_r = u.pow(self.r.clone().get_value());

        // Check 1: u^r = w
        let check1 = u_r == *w;

        // Check 2: u = v^t
        let check2 = *u == v.pow(t_s.get_value());

        check1 && check2
    }
}

pub fn knowledge_of_exponent_protocol<F: Field>(
    prover: &Prover<F>,
    verifier: &Verifier<F>,
) -> bool {
    // Step 1 & 2: Verifier generates a challenge
    let max_degree = std::cmp::max(prover.p.degree(), prover.h.degree()) as usize;
    let (alpha_powers, alpha_prime_powers) = verifier.generate_challenge(max_degree + 1);

    // Step 3: Prover computes and sends u = g^p, v = g^h, and w = g^p'
    let (u, v, w) = prover.compute_values(&alpha_powers, &alpha_prime_powers);

    // Step 4 & 5: Verifier checks whether u^r = w and u = v^t
    verifier.verify(&u, &v, &w)
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::modules::field::{FiniteFieldElement, ModEIP197};
    use crate::modules::polynomial::Polynomial;
    use crate::modules::ring::Ring;

    #[test]
    fn test_kea_protocol() {
        const GENERATOR: i128 = 5; // A primitive root modulo MODULUS

        type F = FiniteFieldElement<ModEIP197>;

        // Create polynomials P(x) and T(x)
        let p = Polynomial::from_monomials(&[F::from_value(1), F::from_value(2), F::from_value(3)]);
        let t = Polynomial::from_monomials(&[F::from_value(1), F::from_value(2)]);

        // Honest protocol
        let honest_prover = Prover::new(p.clone(), t.clone());
        let verifier = Verifier::new(t.clone(), GENERATOR);

        let honest_result = knowledge_of_exponent_protocol(&honest_prover, &verifier);
        println!("Honest protocol result: {}", honest_result);
    }
}
