use std::collections::HashMap;

use crate::modules::field::Field;
use crate::modules::polynomial::Polynomial;

pub struct Prover<F: Field> {
    pub p: Polynomial<F>,
    pub t: Polynomial<F>,
    pub h: Polynomial<F>,
}

pub struct Verifier<F: Field> {
    pub t: Polynomial<F>,
    pub known_roots: Vec<F>,
}

impl<F: Field> Prover<F> {
    pub fn new(p: Polynomial<F>, t: Polynomial<F>) -> Self {
        let h = p.clone() / t.clone();
        Prover { p, t, h }
    }

    pub fn compute_all_values(&self, modulus: i128) -> (HashMap<F, F>, HashMap<F, F>) {
        let mut h_values = HashMap::new();
        let mut p_values = HashMap::new();

        for i in 0..modulus {
            let x = F::from_value(i);
            h_values.insert(x.clone(), self.h.eval(&x));
            p_values.insert(x.clone(), self.p.eval(&x));
        }

        (h_values, p_values)
    }
}

impl<F: Field> Verifier<F> {
    pub fn new(known_roots: Vec<F>) -> Self {
        let t = Polynomial::from_monomials(&known_roots);
        Verifier { t, known_roots }
    }

    pub fn verify(&self, h_values: &HashMap<F, F>, p_values: &HashMap<F, F>) -> bool {
        for (x, h_x) in h_values {
            let t_x = self.t.eval(x);
            let p_x = p_values.get(x).unwrap();
            if h_x.clone() * t_x != *p_x {
                return false;
            }
        }
        true
    }
}

pub fn naive_protocol<F: Field>(prover: &Prover<F>, verifier: &Verifier<F>, modulus: i128) -> bool {
    // Step 1: Verifier sends all possible values (implicitly done by Prover computing all values)

    // Step 2: Prover computes and sends all possible outputs
    let (h_values, p_values) = prover.compute_all_values(modulus);

    // Step 3: Verifier checks whether H(a)T(a) = P(a) holds for any a in F
    verifier.verify(&h_values, &p_values)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        define_myzkp_modulus_type,
        modules::field::{FiniteFieldElement, ModulusValue},
        modules::polynomial::Polynomial,
    };
    use num_bigint::BigInt;
    use std::str::FromStr;

    #[test]
    fn test_naive_protocol() {
        define_myzkp_modulus_type!(Mod31, "31");

        type F = FiniteFieldElement<Mod31>;

        // Create a polynomial P(x) = (x - 1)(x - 2)(x - 3)(x - 4)(x - 5)
        let roots: Vec<F> = (1..=5).map(|i| F::from_value(i)).collect();
        let p = Polynomial::from_monomials(&roots);

        // Verifier knows the first 3 roots
        let known_roots: Vec<F> = roots[0..3].to_vec();
        let t = Polynomial::from_monomials(&known_roots);

        let prover = Prover::new(p.clone(), t.clone());
        let verifier = Verifier::new(known_roots);

        let result = naive_protocol(&prover, &verifier, 31);
        assert!(result);
    }
}
