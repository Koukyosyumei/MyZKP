use std::ops::Sub;

use num_traits::{One, Zero};

use crate::modules::algebra::curve::bn128::{optimal_ate_pairing, FqOrder, G1Point, G2Point};
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::ring::Ring;

pub struct PublicKeyKZG {
    pub powers_1: Vec<G1Point>,
    pub powers_2: Vec<G2Point>,
}

pub type CommitmentKZG = G1Point;

pub struct ProofKZG {
    pub y: FqOrder,
    pub w: G1Point,
}

pub struct BatchProofKZG {
    pub ys: Vec<FqOrder>,
    pub w: G1Point,
}

pub type ProofDegreeBound = G1Point;

pub fn setup_kzg(g1: &G1Point, g2: &G2Point, n: usize) -> PublicKeyKZG {
    let alpha = FqOrder::random_element(&[]);

    let mut powers_1 = Vec::with_capacity(n);
    let mut alpha_power = FqOrder::one();
    for _ in 0..1 + n {
        powers_1.push(g1.mul_ref(alpha_power.clone().get_value()));
        alpha_power = alpha_power * alpha.clone();
    }

    let powers_2 = vec![g2.clone(), g2.mul_ref(alpha.get_value())];

    PublicKeyKZG { powers_1, powers_2 }
}

pub fn setup_kzg_with_full_g2(g1: &G1Point, g2: &G2Point, n: usize) -> PublicKeyKZG {
    let alpha = FqOrder::random_element(&[]);

    let mut powers_1 = Vec::with_capacity(n);
    let mut powers_2 = Vec::with_capacity(n);
    let mut alpha_power = FqOrder::one();
    for _ in 0..1 + n {
        powers_1.push(g1.mul_ref(alpha_power.clone().get_value()));
        powers_2.push(g2.mul_ref(alpha_power.clone().get_value()));
        alpha_power = alpha_power * alpha.clone();
    }

    PublicKeyKZG { powers_1, powers_2 }
}

pub fn commit_kzg(p: &Polynomial<FqOrder>, pk: &PublicKeyKZG) -> CommitmentKZG {
    p.eval_with_powers_on_curve(&pk.powers_1)
}

pub fn open_kzg(p: &Polynomial<FqOrder>, z: &FqOrder, pk: &PublicKeyKZG) -> ProofKZG {
    let y = p.eval(z);
    let y_poly = Polynomial {
        coef: (&[y.clone()]).to_vec(),
    };
    let q = (p - &y_poly) / Polynomial::<FqOrder>::from_monomials(&[z.clone()]);

    ProofKZG {
        y: y,
        w: q.eval_with_powers_on_curve(&pk.powers_1),
    }
}

pub fn batch_open_kzg(
    p: &Polynomial<FqOrder>,
    zs: &Vec<FqOrder>,
    pk: &PublicKeyKZG,
) -> BatchProofKZG {
    let ys = zs.iter().map(|z| p.eval(z)).collect::<Vec<_>>();
    let ip = Polynomial::interpolate(zs, &ys);
    let denominator = Polynomial::from_monomials(zs);
    let q = (p - &ip) / denominator;

    BatchProofKZG {
        ys: ys,
        w: q.eval_with_powers_on_curve(&pk.powers_1),
    }
}

pub fn verify_kzg(z: &FqOrder, c: &CommitmentKZG, proof: &ProofKZG, pk: &PublicKeyKZG) -> bool {
    let g1 = &pk.powers_1[0];
    let g2 = &pk.powers_2[0];
    let g2_s = &pk.powers_2[1];
    let g2_z = g2.mul_ref(z.clone().get_value());
    let g2_s_minus_z = g2_s.clone() - g2_z;

    let e1 = optimal_ate_pairing(&proof.w, &g2_s_minus_z);
    let e2 = optimal_ate_pairing(&g1, &g2);
    let e3 = optimal_ate_pairing(&c, &g2);

    e3 == e1 * (e2.pow(proof.y.get_value()))
}

pub fn batch_verify_kzg(
    zs: &Vec<FqOrder>,
    c: &CommitmentKZG,
    proof: &BatchProofKZG,
    pk: &PublicKeyKZG,
) -> bool {
    let g2 = &pk.powers_2[0];
    let ip = Polynomial::interpolate(zs, &proof.ys);
    let denominator = Polynomial::from_monomials(zs);
    let g1_ip = ip.eval_with_powers_on_curve(&pk.powers_1);
    let g2_z = denominator.eval_with_powers_on_curve(&pk.powers_2);

    let e1 = optimal_ate_pairing(&proof.w, &g2_z);
    let e2 = optimal_ate_pairing(&(c.clone() - g1_ip), g2);
    e1 == e2
}

pub fn prove_degree_bound(
    p: &Polynomial<FqOrder>,
    pk: &PublicKeyKZG,
    d: usize,
) -> ProofDegreeBound {
    let max_d = pk.powers_1.len() - 1;
    let mut q_coef = (0..(max_d + 1 - d))
        .map(|_| FqOrder::zero())
        .collect::<Vec<_>>();
    q_coef[max_d - d] = FqOrder::one();
    let q = Polynomial { coef: q_coef };
    let r = p * &q;
    r.eval_with_powers_on_curve(&pk.powers_1)
}

pub fn verify_degree_bound(
    c: &CommitmentKZG,
    proof: &ProofDegreeBound,
    pk: &PublicKeyKZG,
    d: usize,
) -> bool {
    let max_d = pk.powers_1.len() - 1;
    optimal_ate_pairing(proof, &pk.powers_2[0]) == optimal_ate_pairing(c, &pk.powers_2[max_d - d])
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::modules::algebra::curve::bn128::BN128;

    #[test]
    fn test_kzg() {
        let g1 = BN128::generator_g1();
        let g2 = BN128::generator_g2();

        let p = Polynomial::from_monomials(&[
            FqOrder::from_value(-1),
            FqOrder::from_value(-2),
            FqOrder::from_value(-3),
        ]);

        let pk = setup_kzg(&g1, &g2, 3);

        // Commitment
        let c = commit_kzg(&p, &pk);

        // Proof
        let z = FqOrder::from_value(5);
        let proof = open_kzg(&p, &z, &pk);

        // Verification
        let flag = verify_kzg(&z, &c, &proof, &pk);
        assert!(flag);
    }

    #[test]
    fn test_batch_kzg() {
        let g1 = BN128::generator_g1();
        let g2 = BN128::generator_g2();

        let p = Polynomial::from_monomials(&[
            FqOrder::from_value(-1),
            FqOrder::from_value(-2),
            FqOrder::from_value(-3),
        ]);

        let pk = setup_kzg_with_full_g2(&g1, &g2, 3);

        // Commitment
        let c = commit_kzg(&p, &pk);

        // Proof
        let zs = vec![FqOrder::from_value(5), FqOrder::from_value(7)];
        let mut proof = batch_open_kzg(&p, &zs, &pk);

        // Verification
        let flag = batch_verify_kzg(&zs, &c, &proof, &pk);
        assert!(flag);

        // Verification
        proof.ys[0] += FqOrder::one();
        let flag = batch_verify_kzg(&zs, &c, &proof, &pk);
        assert!(!flag);
    }

    #[test]
    fn test_degree_bound() {
        let g1 = BN128::generator_g1();
        let g2 = BN128::generator_g2();

        let p = Polynomial::from_monomials(&[
            FqOrder::from_value(-1),
            FqOrder::from_value(-2),
            FqOrder::from_value(-3),
        ]);

        let pk = setup_kzg_with_full_g2(&g1, &g2, 4);

        // Commitment
        let c = commit_kzg(&p, &pk);

        let proof = prove_degree_bound(&p, &pk, 3);
        assert!(verify_degree_bound(&c, &proof, &pk, 3));
    }
}
