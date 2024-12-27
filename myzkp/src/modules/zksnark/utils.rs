use num_traits::Zero;

use crate::modules::algebra::curve::curve::{EllipticCurve, EllipticCurvePoint};
use crate::modules::algebra::field::Field;
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::arithmetization::qap::QAP;

pub fn generate_challenge_vec<F1: Field, F2: Field, E: EllipticCurve>(
    point: &EllipticCurvePoint<F1, E>,
    poly_vec: &[Polynomial<F2>],
    s: &F2,
) -> Vec<EllipticCurvePoint<F1, E>> {
    poly_vec
        .iter()
        .map(|poly| point * poly.eval(s).sanitize().get_value())
        .collect()
}

pub fn generate_alpha_challenge_vec<F1: Field, F2: Field, E: EllipticCurve>(
    point: &EllipticCurvePoint<F1, E>,
    poly_vec: &[Polynomial<F2>],
    s: &F2,
    alpha: &F2,
) -> Vec<EllipticCurvePoint<F1, E>> {
    poly_vec
        .iter()
        .map(|poly| point * (alpha.mul_ref(&poly.eval(s).sanitize())).get_value())
        .collect()
}

pub fn generate_s_powers<F1: Field, F2: Field, E: EllipticCurve>(
    point: &EllipticCurvePoint<F1, E>,
    s: &F2,
    m: usize,
) -> Vec<EllipticCurvePoint<F1, E>> {
    let mut powers = Vec::with_capacity(m + 1);
    let mut current = F2::one();
    for _ in 0..=m {
        powers.push(point * current.get_value());
        current = current * s.clone();
    }
    powers
}

pub fn accumulate_curve_points<F1: Field, F2: Field, E: EllipticCurve>(
    g_vec: &[EllipticCurvePoint<F1, E>],
    assignment: &[F2],
) -> EllipticCurvePoint<F1, E> {
    g_vec.iter().zip(assignment.iter()).fold(
        EllipticCurvePoint::<F1, E>::point_at_infinity(),
        |acc, (g, &ref a)| acc + g * a.get_value(),
    )
}

pub fn accumulate_polynomials<F: Field>(
    poly_vec: &[Polynomial<F>],
    assignment: &[F],
) -> Polynomial<F> {
    poly_vec
        .iter()
        .zip(assignment.iter())
        .fold(Polynomial::<F>::zero(), |acc, (poly, &ref a)| {
            acc + poly.clone() * a.clone()
        })
}

pub fn get_h<F: Field>(qap: &QAP<F>, assignment: &Vec<F>) -> Polynomial<F> {
    let ell = accumulate_polynomials(&qap.ell_i_vec, assignment);
    let r = accumulate_polynomials(&qap.r_i_vec, assignment);
    let o = accumulate_polynomials(&qap.o_i_vec, assignment);
    (ell * r - o) / qap.t.clone()
}
