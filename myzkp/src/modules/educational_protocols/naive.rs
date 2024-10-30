use std::collections::HashMap;

use crate::modules::field::Field;
use crate::modules::polynomial::Polynomial;

struct Prover<F: Field> {
    p: Polynomial<F>,
    t: Polynomial<F>,
    h: Polynomial<F>,
}

struct Verifier<F: Field> {
    t: Polynomial<F>,
    known_roots: Vec<F>,
}

impl<F: Field> Prover<F> {
    fn new(p: Polynomial<F>, t: Polynomial<F>) -> Self {
        let h = p.clone() / t.clone();
        Prover { p, t, h }
    }

    fn compute_all_values(&self, modulus: i128) -> (HashMap<F, F>, HashMap<F, F>) {
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
    fn new(known_roots: Vec<F>) -> Self {
        let t = Polynomial::from_monomials(&known_roots);
        Verifier { t, known_roots }
    }

    fn verify(&self, h_values: &HashMap<F, F>, p_values: &HashMap<F, F>) -> bool {
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

fn naive_protocol<F: Field>(prover: &Prover<F>, verifier: &Verifier<F>, modulus: i128) -> bool {
    // Step 1: Verifier sends all possible values (implicitly done by Prover computing all values)

    // Step 2: Prover computes and sends all possible outputs
    let (h_values, p_values) = prover.compute_all_values(modulus);

    // Step 3: Verifier checks whether H(a)T(a) = P(a) holds for any a in F
    verifier.verify(&h_values, &p_values)
}
