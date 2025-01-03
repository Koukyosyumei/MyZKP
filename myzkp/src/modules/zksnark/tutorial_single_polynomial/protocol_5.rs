use num_bigint::ToBigInt;

use crate::modules::algebra::field::Field;
use crate::modules::algebra::polynomial::Polynomial;

pub struct Prover5<F: Field> {
    pub p: Polynomial<F>,
    pub t: Polynomial<F>,
    pub h: Polynomial<F>,
}

pub struct Verifier5<F: Field> {
    t: Polynomial<F>,
    s: F,
    r: F,
    g: F,
}

impl<F: Field> Prover5<F> {
    pub fn new(p: Polynomial<F>, t: Polynomial<F>) -> Self {
        let h = p.clone() / t.clone();
        Prover5 { p, t, h }
    }

    pub fn compute_values(&self, s_powers: &[F], s_prime_powers: &[F]) -> (F, F, F) {
        let delta = F::random_element(&[]);

        let g_p = self.p.eval_with_powers(s_powers);
        let g_h = self.h.eval_with_powers(s_powers);
        let g_p_prime = self.p.eval_with_powers(s_prime_powers);

        let u_prime = g_p.pow(delta.get_value());
        let v_prime = g_h.pow(delta.get_value());
        let w_prime = g_p_prime.pow(delta.get_value());

        (u_prime, v_prime, w_prime)
    }
}

impl<F: Field> Verifier5<F> {
    pub fn new(t: Polynomial<F>, generator: i128) -> Self {
        let s = F::random_element(&[]);
        let r = F::random_element(&[]);
        let g = F::from_value(generator);
        Verifier5 { t, s, r, g }
    }

    pub fn generate_challenge(&self, max_degree: usize) -> (Vec<F>, Vec<F>) {
        let mut s_powers = vec![];
        let mut s_prime_powers = vec![];

        for i in 0..(max_degree + 1) {
            s_powers.push(
                self.g
                    .pow(self.s.clone().pow_m1(i.to_bigint().unwrap()).get_value()),
            );
            s_prime_powers.push(s_powers.last().unwrap().pow(self.r.get_value()));
        }

        (s_powers, s_prime_powers)
    }

    pub fn verify(&self, u_prime: &F, v_prime: &F, w_prime: &F) -> bool {
        let t_s = self.t.eval_m1(&self.s);
        let u_prime_r = u_prime.pow(self.r.clone().get_value());

        // Check 1: u'^r = w'
        let check1 = u_prime_r == *w_prime;

        // Check 2: u' = v'^t
        let check2 = *u_prime == v_prime.pow(t_s.get_value());

        check1 && check2
    }
}

pub fn zk_protocol<F: Field>(prover: &Prover5<F>, verifier: &Verifier5<F>) -> bool {
    // Step 1 & 2: Verifier5 generates a challenge
    let max_degree = std::cmp::max(prover.p.degree(), prover.h.degree()) as usize;
    let (s_powers, s_prime_powers) = verifier.generate_challenge(max_degree + 1);

    // Step 3 & 4: Prover5 computes and sends u' = (g^p)^δ, v' = (g^h)^δ, and w' = (g^p')^δ
    let (u_prime, v_prime, w_prime) = prover.compute_values(&s_powers, &s_prime_powers);

    // Step 5 & 6: Verifier5 checks whether u'^r = w' and u' = v'^t
    verifier.verify(&u_prime, &v_prime, &w_prime)
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::modules::algebra::field::{FiniteFieldElement, ModEIP197};
    use crate::modules::algebra::polynomial::Polynomial;
    use crate::modules::algebra::ring::Ring;

    #[test]
    fn test_kea_protocol() {
        const GENERATOR: i128 = 5; // A primitive root modulo MODULUS

        type F = FiniteFieldElement<ModEIP197>;

        // Create polynomials P(x) and T(x)
        let p = Polynomial::from_monomials(&[F::from_value(1), F::from_value(2), F::from_value(3)]);
        let t = Polynomial::from_monomials(&[F::from_value(1), F::from_value(2)]);

        // ZK protocol
        let prover = Prover5::new(p.clone(), t.clone());
        let verifier = Verifier5::new(t.clone(), GENERATOR);

        let result = zk_protocol(&prover, &verifier);
        println!("ZK protocol result: {}", result);
    }
}
