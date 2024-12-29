use num_traits::Zero;

use crate::modules::algebra::curve::curve::{EllipticCurve, EllipticCurvePoint};
use crate::modules::algebra::field::Field;
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::arithmetization::qap::QAP;

/// Generates a vector of elliptic curve points by evaluating each polynomial in the input vector
/// at a specified point `s` and scaling the input curve point by the evaluation result.
///
/// # Arguments
/// * `point` - Reference to the base elliptic curve point.
/// * `poly_vec` - A slice of polynomials to be evaluated.
/// * `s` - Reference to the field element at which the polynomials are evaluated.
///
/// # Returns
/// A vector of elliptic curve points obtained by scaling the input `point` by the evaluated results.
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

/// Generates a vector of elliptic curve points by evaluating polynomials at `s`,
/// multiplying the results by a challenge `alpha`, and scaling the input curve point accordingly.
///
/// # Arguments
/// * `point` - Reference to the base elliptic curve point.
/// * `poly_vec` - A slice of polynomials to be evaluated.
/// * `s` - Reference to the field element at which the polynomials are evaluated.
/// * `alpha` - Reference to the challenge field element used for scaling.
///
/// # Returns
/// A vector of elliptic curve points scaled by the evaluated and challenged results.
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

/// Computes powers of a field element `s` and scales the given elliptic curve point accordingly.
///
/// # Arguments
/// * `point` - Reference to the base elliptic curve point.
/// * `s` - Reference to the field element used to compute powers.
/// * `m` - The maximum power of `s` to compute.
///
/// # Returns
/// A vector of elliptic curve points corresponding to `point` scaled by powers of `s`.
pub fn generate_s_powers<F1: Field, F2: Field, E: EllipticCurve>(
    point: &EllipticCurvePoint<F1, E>,
    s: &F2,
    m: usize,
) -> Vec<EllipticCurvePoint<F1, E>> {
    let mut powers = Vec::with_capacity(m + 1);
    let mut current = F2::one();
    for _ in 0..=m {
        powers.push(point * current.get_value());
        current *= s.clone();
    }
    powers
}

/// Accumulates a weighted sum of elliptic curve points based on the given assignments.
///
/// # Arguments
/// * `g_vec` - A slice of elliptic curve points.
/// * `assignment` - A slice of field elements serving as weights.
///
/// # Returns
/// The accumulated elliptic curve point representing the weighted sum.
pub fn accumulate_curve_points<F1: Field, F2: Field, E: EllipticCurve>(
    g_vec: &[EllipticCurvePoint<F1, E>],
    assignment: &[F2],
) -> EllipticCurvePoint<F1, E> {
    g_vec.iter().zip(assignment.iter()).fold(
        EllipticCurvePoint::<F1, E>::point_at_infinity(),
        |acc, (g, &ref a)| acc + g * a.get_value(),
    )
}

/// Computes a weighted sum of polynomials based on the given assignments.
///
/// # Arguments
/// * `poly_vec` - A slice of polynomials.
/// * `assignment` - A slice of field elements serving as weights.
///
/// # Returns
/// A polynomial representing the weighted sum of the input polynomials.
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

/// Computes the quotient polynomial `h` for a Quadratic Arithmetic Program (QAP).
///
/// The quotient `h` is calculated as:
/// \[ h(x) = \frac{\ell(x) \cdot r(x) - o(x)}{t(x)} \]
/// where `\ell(x)`, `r(x)`, and `o(x)` are accumulated polynomials over `qap.ell_i_vec`,
/// `qap.r_i_vec`, and `qap.o_i_vec` respectively.
///
/// # Arguments
/// * `qap` - Reference to a QAP structure containing the polynomials `t`, `ell_i_vec`, `r_i_vec`, and `o_i_vec`.
/// * `assignment` - A vector of field elements used for accumulating polynomials.
///
/// # Returns
/// The resulting quotient polynomial `h`.
pub fn get_h<F: Field>(qap: &QAP<F>, assignment: &Vec<F>) -> Polynomial<F> {
    let ell = accumulate_polynomials(&qap.ell_i_vec, assignment);
    let r = accumulate_polynomials(&qap.r_i_vec, assignment);
    let o = accumulate_polynomials(&qap.o_i_vec, assignment);
    (ell * r - o) / qap.t.clone()
}
