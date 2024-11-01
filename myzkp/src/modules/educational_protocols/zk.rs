use crate::modules::field::Field;
use crate::modules::polynomial::Polynomial;
use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_traits::Zero;
use std::str::FromStr;

struct Prover<F: Field> {
    p: Polynomial<F>,
    t: Polynomial<F>,
    h: Polynomial<F>,
}

struct Verifier<F: Field> {
    t: Polynomial<F>,
    s: F,
    r: F,
    g: F,
}

impl<F: Field> Prover<F> {
    fn new(p: Polynomial<F>, t: Polynomial<F>) -> Self {
        let h = p.clone() / t.clone();
        Prover { p, t, h }
    }

    fn compute_values(&self, g: &F, alpha_powers: &[F], alpha_prime_powers: &[F]) -> (F, F, F) {
        let mut rng = rand::thread_rng();
        let delta = FiniteFieldElement::<MODULUS>::from(rng.gen_range(0..MODULUS));

        let g_p = self.p.eval_with_powers(g, alpha_powers);
        let g_h = self.h.eval_with_powers(g, alpha_powers);
        let g_p_prime = self.p.eval_with_powers(g, alpha_prime_powers);

        let u_prime = g_p.pow(delta.value);
        let v_prime = g_h.pow(delta.value);
        let w_prime = g_p_prime.pow(delta.value);

        (u_prime, v_prime, w_prime)
    }
}

impl<F: Field> Verifier<F> {
    fn new(t: Polynomial<F>) -> Self {
        let mut rng = rand::thread_rng();
        let s = FiniteFieldElement::<MODULUS>::from(rng.gen_range(0..MODULUS));
        let r = FiniteFieldElement::<MODULUS>::from(rng.gen_range(0..MODULUS));
        let g = FiniteFieldElement::<MODULUS>::from(GENERATOR);
        Verifier { t, s, r, g }
    }

    fn generate_challenge(&self, max_degree: usize) -> (F, Vec<F>, Vec<F>) {
        let alpha = self.g.pow(self.s.clone().value);
        let alpha_prime = alpha.pow(self.r.clone().value);

        let mut alpha_powers = vec![alpha.clone()];
        let mut alpha_prime_powers = vec![alpha_prime.clone()];

        for _ in 1..max_degree {
            alpha_powers.push(alpha_powers.last().unwrap() * alpha.clone());
            alpha_prime_powers.push(alpha_prime_powers.last().unwrap() * alpha_prime.clone());
        }

        (self.g.clone(), alpha_powers, alpha_prime_powers)
    }

    fn verify(&self, u_prime: &F, v_prime: &F, w_prime: &F) -> bool {
        let t_s = self.t.eval(&self.s);
        let u_prime_r = u_prime.pow(self.r.clone().value);

        // Check 1: u'^r = w'
        let check1 = u_prime_r == *w_prime;

        // Check 2: u' = v'^t
        let check2 = *u_prime == v_prime.pow(t_s.value);

        check1 && check2
    }
}

fn zk_snark_protocol<F: Field>(prover: &Prover<F>, verifier: &Verifier<F>) -> bool {
    // Step 1 & 2: Verifier generates a challenge
    let max_degree = std::cmp::max(prover.p.degree(), prover.h.degree()) as usize;
    let (g, alpha_powers, alpha_prime_powers) = verifier.generate_challenge(max_degree + 1);

    // Step 3 & 4: Prover computes and sends u' = (g^p)^δ, v' = (g^h)^δ, and w' = (g^p')^δ
    let (u_prime, v_prime, w_prime) = prover.compute_values(&g, &alpha_powers, &alpha_prime_powers);

    // Step 5 & 6: Verifier checks whether u'^r = w' and u' = v'^t
    verifier.verify(&u_prime, &v_prime, &w_prime)
}

fn main() {
    type F = FiniteFieldElement<MODULUS>;

    // Create polynomials P(x) and T(x)
    let p = Polynomial::from_monomials(&[F::from(1), F::from(2), F::from(3)]);
    let t = Polynomial::from_monomials(&[F::from(1), F::from(2)]);

    // ZK-SNARK protocol
    let prover = Prover::new(p.clone(), t.clone());
    let verifier = Verifier::new(t.clone());

    let result = zk_snark_protocol(&prover, &verifier);
    println!("ZK-SNARK protocol result: {}", result);
}
