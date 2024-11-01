use crate::modules::field::{Field, FiniteFieldElement};
use crate::modules::polynomial::Polynomial;
use rand::Rng;

impl<F: Field> Polynomial<F> {
    fn eval_with_powers(&self, g: &F, alpha_powers: &[F]) -> F {
        let mut result = F::one();
        for (i, coef) in self.poly.iter().enumerate() {
            if i == 0 {
                result = result * g.pow(coef.clone().get_value());
            } else {
                result = result * alpha_powers[i - 1].pow(coef.clone().get_value());
            }
        }
        result
    }
}

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
        let g_p = self.p.eval_with_powers(g, alpha_powers);
        let g_h = self.h.eval_with_powers(g, alpha_powers);
        let g_p_prime = self.p.eval_with_powers(g, alpha_prime_powers);
        (g_p, g_h, g_p_prime)
    }
}

impl<F: Field> Verifier<F> {
    fn new(t: Polynomial<F>, generatir: i128) -> Self {
        let mut rng = rand::thread_rng();
        let s = F::from_value(rng.gen_bigint_range(
            &BigInt::zero(),
            &BigInt::from_str("4835703278458516698824704").unwrap(), // 2^82
        ));
        let r = F::from_value(rng.gen_bigint_range(
            &BigInt::zero(),
            &BigInt::from_str("4835703278458516698824704").unwrap(), // 2^82
        ));
        let g = F::from_value(generator);
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

    fn verify(&self, u: &F, v: &F, w: &F) -> bool {
        let t_s = self.t.eval(&self.s);
        let u_r = u.pow(self.r.clone().value);

        // Check 1: u^r = w
        let check1 = u_r == *w;

        // Check 2: u = v^t
        let check2 = *u == v.pow(t_s.value);

        check1 && check2
    }
}

fn knowledge_of_exponent_protocol<F: Field>(prover: &Prover<F>, verifier: &Verifier<F>) -> bool {
    // Step 1 & 2: Verifier generates a challenge
    let max_degree = std::cmp::max(prover.p.degree(), prover.h.degree()) as usize;
    let (g, alpha_powers, alpha_prime_powers) = verifier.generate_challenge(max_degree + 1);

    // Step 3: Prover computes and sends u = g^p, v = g^h, and w = g^p'
    let (u, v, w) = prover.compute_values(&g, &alpha_powers, &alpha_prime_powers);

    // Step 4 & 5: Verifier checks whether u^r = w and u = v^t
    verifier.verify(&u, &v, &w)
}

fn main() {
    type F = FiniteFieldElement<MODULUS>;

    // Create polynomials P(x) and T(x)
    let p = Polynomial::from_monomials(&[F::from(1), F::from(2), F::from(3)]);
    let t = Polynomial::from_monomials(&[F::from(1), F::from(2)]);

    // Honest protocol
    let honest_prover = Prover::new(p.clone(), t.clone());
    let verifier = Verifier::new(t.clone());

    let honest_result = knowledge_of_exponent_protocol(&honest_prover, &verifier);
    println!("Honest protocol result: {}", honest_result);
}
