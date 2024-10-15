use crate::modules::field::ModuloFieldElement;
use crate::modules::polynomial::Polynomial;

pub struct SNARK {
    g: ModuloFieldElement,
    t: Polynomial,
    n: i32,
}

pub struct ProofKey {
    alpha: ModuloFieldElement,
    alpha_seq: Vec<ModuloFieldElement>,
    alpha_prime: ModuloFieldElement,
    alpha_prime_seq: Vec<ModuloFieldElement>,
}

pub struct VerificationKey {
    gr: ModuloFieldElement,
    gts: ModuloFieldElement,
}

pub struct Proof {
    u_prime: ModuloFieldElement,
    v_prime: ModuloFieldElement,
    w_prime: ModuloFieldElement,
}

/// Function simulating bilinear pairing on two group elements.
//fn simulate_bilinear_pairing(g1: ModuloFieldElement, g2: ModuloFieldElement) -> ModuloFieldElement {
//    // Pairing function e(g1, g2) = g1.element^g2.element mod modulus
//    g1.element.pow(g2.element.value.to_u64_digits()[0]) // Simplified pairing logic
//}

impl SNARK {
    pub fn trusted_setup(&self) -> (ProofKey, VerificationKey) {
        let excluded_elements: Vec<ModuloFieldElement> = vec![];
        let s = ModuloFieldElement::random_element(&excluded_elements);
        let r = ModuloFieldElement::random_element(&excluded_elements);

        let alpha = self.g.pow(s.clone());
        let alpha_prime = self.g.pow(s.clone() * r.clone()); // g^(sr)

        let proof_key = ProofKey {
            alpha,
            alpha_seq: (1..self.n)
                .map(|i| self.g.pow(s.mul_scalar(i.into())))
                .collect(),
            alpha_prime,
            alpha_prime_seq: (1..self.n)
                .map(|i| self.g.pow((s.clone() * r.clone()).mul_scalar(i.into())))
                .collect(),
        };

        let verifier_key = VerificationKey {
            gr: self.g.pow(r),                // g^r
            gts: self.g.pow(self.t.eval(&s)), // g^T(s)
        };

        (proof_key, verifier_key)
    }

    pub fn prove(&self, proof_key: &ProofKey, p: &Polynomial, h: &Polynomial) -> Proof {
        let excluded_elements: Vec<ModuloFieldElement> = vec![];
        let delta = ModuloFieldElement::random_element(&excluded_elements);

        let mut gp = self.g.pow(p.poly[0].clone());
        for i in 1..p.poly.len() {
            gp = gp * (proof_key.alpha_seq[i - 1]).pow(p.poly[i].clone());
        }
        let u_prime = gp.pow(delta.clone());

        let mut gh = self.g.pow(h.poly[0].clone());
        for i in 1..h.poly.len() {
            gh = gh * (proof_key.alpha_seq[i - 1]).pow(h.poly[i].clone());
        }
        let v_prime = gh.pow(delta.clone());

        let mut gp_prime = self.g.pow(p.poly[0].clone());
        for i in 1..p.poly.len() {
            gp_prime = gp_prime * (proof_key.alpha_prime_seq[i - 1]).pow(p.poly[i].clone());
        }
        let w_prime = gp_prime.pow(delta.clone());

        Proof {
            u_prime: u_prime,
            v_prime: v_prime,
            w_prime: w_prime,
        }
    }

    //pub fn verify(&self, verification_key: &VerificationKey, proof: Proof) -> bool {}
}
