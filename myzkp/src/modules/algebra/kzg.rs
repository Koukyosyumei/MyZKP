use num_traits::One;

use crate::modules::algebra::curve::bn128::{optimal_ate_pairing, FqOrder, G1Point, G2Point};
use crate::modules::algebra::field::Field;
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::ring::Ring;

pub struct PublicKey {
    alpha_1: Vec<G1Point>,
    alpha_2: Vec<G2Point>,
}

pub fn setup_kzg(g1: &G1Point, g2: &G2Point, n: usize) -> PublicKey {
    let s = FqOrder::random_element(&[]);

    let mut alpha_1 = Vec::with_capacity(n);
    let mut alpha_2 = Vec::with_capacity(n);

    let mut s_power = FqOrder::one();
    for _ in 0..1 + n {
        alpha_1.push(g1.mul_ref(s_power.clone().get_value()));
        alpha_2.push(g2.mul_ref(s_power.clone().get_value()));
        s_power = s_power * s.clone();
    }

    PublicKey { alpha_1, alpha_2 }
}

pub fn commit_kzg(p: &Polynomial<FqOrder>, pk: &PublicKey) -> G1Point {
    p.eval_with_powers_on_curve(&pk.alpha_1)
}

pub fn open_kzg(p: &Polynomial<FqOrder>, z: &FqOrder, pk: &PublicKey) -> G1Point {
    let y = p.eval(z);
    let y_poly = Polynomial {
        coef: (&[y]).to_vec(),
    };

    let q = (p - &y_poly) / Polynomial::<FqOrder>::from_monomials(&[z.clone()]);
    q.eval_with_powers_on_curve(&pk.alpha_1)
}

pub fn verify_kzg(z: &FqOrder, y: &FqOrder, c: &G1Point, w: &G1Point, pk: &PublicKey) -> bool {
    let g1 = &pk.alpha_1[0];
    let g2 = &pk.alpha_2[0];
    let g2_s = &pk.alpha_2[1];
    let g2_z = g2.mul_ref(z.clone().get_value());
    let g2_s_minus_z = g2_s.clone() - g2_z;

    let e1 = optimal_ate_pairing(&w, &g2_s_minus_z);
    let e2 = optimal_ate_pairing(&g1, &g2);
    let e3 = optimal_ate_pairing(&c, &g2);

    e3 == e1 * (e2.pow(y.get_value()))
}
