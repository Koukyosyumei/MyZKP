use crate::modules::field::Field;
use crate::modules::polynomial::Polynomial;
use num_bigint::ToBigInt;

impl<F: Field> Polynomial<F> {
    fn eval_with_powers(&self, g: &F, alpha_powers: &[F]) -> F {
        let mut result = F::one();
        println!("f: {}", result);
        for (i, coef) in self.poly.iter().enumerate() {
            if i == 0 {
                result = result * g.pow(coef.clone().get_value());
                println!("coef: {}, result: {}", coef.get_value(), result);
            } else {
                result = result * alpha_powers[i - 1].pow(coef.clone().get_value());
                println!("coef: {}, result: {}", coef.get_value(), result);
            }
        }
        result
    }
}

pub struct Prover<F: Field> {
    p: Polynomial<F>,
    t: Polynomial<F>,
    h: Polynomial<F>,
}

pub struct Verifier<F: Field> {
    t: Polynomial<F>,
    s: F,
    g: F,
}

impl<F: Field> Prover<F> {
    pub fn new(p: Polynomial<F>, t: Polynomial<F>) -> Self {
        println!("p={}, t={}", p, t);
        let h = p.clone() / t.clone();
        println!("h={}", h);
        Prover { p, t, h }
    }

    pub fn compute_values(&self, g: &F, alpha_powers: &[F]) -> (F, F) {
        let g_p = self.p.eval_with_powers(g, alpha_powers);
        let g_h = self.h.eval_with_powers(g, alpha_powers);
        (g_p, g_h)
    }
}

impl<F: Field> Verifier<F> {
    pub fn new(t: Polynomial<F>, generator: i128) -> Self {
        let s = F::from_value(1477); // 1477 bad, 1476 ok
        let g = F::from_value(generator);
        println!("inv of g: {} ({})", g.inverse(), g);
        Verifier { t, s, g }
    }

    pub fn generate_challenge(&self, max_degree: usize) -> (F, Vec<F>) {
        //let alpha = self.g.pow(self.s.clone().get_value());
        let mut alpha_powers = vec![];
        for i in 1..(max_degree + 1) {
            alpha_powers.push(
                self.g
                    .pow(self.s.clone().pow(i.to_bigint().unwrap()).get_value()),
            );
            println!(
                "aa: {}",
                self.g.pow(
                    (self.s.clone() - F::one() - F::one())
                        .clone()
                        .pow(i.to_bigint().unwrap())
                        .get_value()
                )
            );
            println!(
                "ab: {}",
                self.g.pow(
                    (self.s.clone() - F::one())
                        .clone()
                        .pow(i.to_bigint().unwrap())
                        .get_value()
                )
            );
            println!("ap: {}", alpha_powers.last().unwrap());
            //alpha_powers.push(alpha_powers.last().unwrap().clone() * alpha.clone());
        }
        (self.g.clone(), alpha_powers)
    }

    pub fn verify(&self, u: &F, v: &F) -> bool {
        let t_s = self.t.eval(&self.s);
        println!("s={}, t={}, t_s={}", self.s, self.t, t_s);
        println!("u={}, v.pow.t_s={}", u, &v.pow(t_s.get_value()));
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

    pub fn compute_malicious_values(&self, g: &F, alpha_powers: &[F]) -> (F, F) {
        let g_t = self.t.eval_with_powers(g, alpha_powers);
        let z = F::random_element(&[]);
        let fake_v = g.pow(z.get_value());
        let fake_u = g_t.pow(z.get_value());
        (fake_u, fake_v)
    }
}

pub fn discrete_log_protocol<F: Field>(prover: &Prover<F>, verifier: &Verifier<F>) -> bool {
    // Step 1 & 2: Verifier generates a challenge
    let max_degree = prover.p.degree();
    println!("md: {}", max_degree);
    let (g, alpha_powers) = verifier.generate_challenge(max_degree as usize);
    println!("g: {}", g);
    for a in &alpha_powers {
        println!("a: {}", a);
    }

    // Step 3: Prover computes and sends u = g^p and v = g^h
    let (u, v) = prover.compute_values(&g, &alpha_powers);
    println!("u: {}", u);
    println!("v: {}", v);

    // Step 4: Verifier checks whether u = v^t
    verifier.verify(&u, &v)
}

pub fn malicious_discrete_log_protocol<F: Field>(
    prover: &MaliciousProver<F>,
    verifier: &Verifier<F>,
) -> bool {
    // Step 1 & 2: Verifier generates a challenge
    let max_degree = prover.t.degree() as usize;
    let (g, alpha_powers) = verifier.generate_challenge(max_degree as usize);

    // Step 3: Malicious Prover computes and sends fake u and v
    let (fake_u, fake_v) = prover.compute_malicious_values(&g, &alpha_powers);

    // Step 4: Verifier checks whether u = v^t (which will pass for the fake values)
    verifier.verify(&fake_u, &fake_v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modules::field::FiniteFieldElement;
    use crate::modules::polynomial::Polynomial;

    #[test]
    fn test_dl_protocol() {
        const MODULUS: i128 = 3 * (1 << 30) + 1;
        const GENERATOR: i128 = 5; // A primitive root modulo MODULUS

        type F = FiniteFieldElement<MODULUS>;

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
        //assert!(malicious_result);
    }
}
