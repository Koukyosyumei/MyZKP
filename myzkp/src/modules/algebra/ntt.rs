use std::ops::Mul;

use num_traits::{One, Zero};

use crate::modules::algebra::field::Field;
use crate::modules::algebra::polynomial::Polynomial;

pub fn ntt<F: Field>(primitive_root: &F, values: &Vec<F>) -> Vec<F> {
    assert!(
        values.len() & (values.len() - 1) == 0,
        "cannot compute ntt of non-power-of-two sequence"
    );
    if values.len() <= 1 {
        return values.to_vec();
    }

    assert!(
        primitive_root.pow(values.len()).is_one(),
        "primitive root must be nth root of unity, where n is len(values)"
    );
    assert!(
        !primitive_root.pow(values.len() / 2).is_one(),
        "primitive root is not primitive nth root of unity, where n is len(values)"
    );

    let half = values.len() / 2;
    let odds = ntt(
        &primitive_root.pow(2),
        &values
            .iter()
            .enumerate()
            .filter(|(index, _)| index % 2 == 1)
            .map(|(_, &ref value)| value.clone())
            .collect(),
    );
    let evens = ntt(
        &primitive_root.pow(2),
        &values
            .iter()
            .enumerate()
            .filter(|(index, _)| index % 2 == 0)
            .map(|(_, &ref value)| value.clone())
            .collect(),
    );

    (0..values.len())
        .map(|i| evens[i % half].add_ref(&((primitive_root.pow(i)) * &odds[i % half])))
        .collect()
}

pub fn intt<F>(primitive_root: &F, values: &Vec<F>) -> Vec<F>
where
    F: Field,
{
    if values.len() == 1 {
        return values.to_vec();
    }

    let ninv = F::from_value(values.len()).inverse();
    let transformed_values = ntt(&primitive_root.inverse(), values);
    transformed_values
        .iter()
        .map(|tv| ninv.mul_ref(tv))
        .collect()
}

pub fn fast_multiply<F>(
    lhs: &Polynomial<F>,
    rhs: &Polynomial<F>,
    primitive_root: &F,
    root_order: usize,
) -> Polynomial<F>
where
    F: Field,
{
    assert!(primitive_root.pow(root_order).is_one());
    assert!(!primitive_root.pow(root_order / 2).is_one());

    if lhs.is_zero() || rhs.is_zero() {
        return Polynomial::<F>::zero();
    }

    let mut root = primitive_root.clone();
    let mut order = root_order;
    let degree = lhs.degree() + rhs.degree();

    if degree < 8 {
        return lhs * rhs;
    }

    while degree < (order / 2).try_into().unwrap() {
        root = root.mul_ref(&root);
        order = order / 2;
    }

    let mut lhs_coefficients = lhs.coef.clone();
    let mut rhs_coefficients = rhs.coef.clone();
    while lhs_coefficients.len() < order {
        lhs_coefficients.push(F::zero());
    }
    while rhs_coefficients.len() < order {
        rhs_coefficients.push(F::zero());
    }

    let lhs_codeword = ntt(&root, &lhs_coefficients);
    let rhs_codeword = ntt(&root, &rhs_coefficients);
    let hadamard_product = lhs_codeword
        .iter()
        .zip(rhs_codeword.iter())
        .map(|(el, r)| el.mul_ref(r))
        .collect();
    let product_coefficients = intt(&root, &hadamard_product);

    Polynomial {
        coef: product_coefficients,
    }
}

pub fn fast_zerofier<F>(domain: &Vec<F>, primitive_root: &F, root_order: usize) -> Polynomial<F>
where
    F: Field,
{
    assert!(primitive_root.pow(root_order).is_one());
    assert!(!primitive_root.pow(root_order / 2).is_one());

    if domain.is_empty() {
        return Polynomial::<F>::zero();
    }

    if domain.len().is_one() {
        return Polynomial {
            coef: vec![F::zero() - &domain[0], F::one()],
        };
    }

    let half = domain.len() / 2;

    let left = fast_zerofier(&domain[..half].to_vec(), primitive_root, root_order);
    let right = fast_zerofier(&domain[half..].to_vec(), primitive_root, root_order);

    fast_multiply(&left, &right, primitive_root, root_order)
}

pub fn fast_evaluate<F>(
    polynomial: &Polynomial<F>,
    domain: &Vec<F>,
    primitive_root: &F,
    root_order: usize,
) -> Vec<F>
where
    F: Field,
{
    assert!(primitive_root.pow(root_order).is_one());
    assert!(!primitive_root.pow(root_order / 2).is_one());

    if domain.len().is_zero() {
        return vec![];
    }

    if domain.len().is_one() {
        return vec![polynomial.eval(&domain[0])];
    }

    let half = domain.len() / 2;

    let left_zerofier = fast_zerofier(&domain[..half].to_vec(), primitive_root, root_order);
    let right_zerofier = fast_zerofier(&domain[half..].to_vec(), primitive_root, root_order);

    let mut left = fast_evaluate(
        &(polynomial % &left_zerofier),
        &domain[..half].to_vec(),
        primitive_root,
        root_order,
    );
    let mut right = fast_evaluate(
        &(polynomial % &right_zerofier),
        &domain[half..].to_vec(),
        primitive_root,
        root_order,
    );

    left.append(&mut right);
    left
}

pub fn fast_interpolate<F>(
    domain: &Vec<F>,
    values: &Vec<F>,
    primitive_root: &F,
    root_order: usize,
) -> Polynomial<F>
where
    F: Field,
{
    assert!(primitive_root.pow(root_order).is_one());
    assert!(!primitive_root.pow(root_order / 2).is_one());
    assert_eq!(domain.len(), values.len());

    if domain.len().is_zero() {
        return Polynomial::zero();
    }

    if domain.len().is_one() {
        return Polynomial {
            coef: vec![values[0].clone()],
        };
    }

    let half = domain.len() / 2;

    let left_zerofier = fast_zerofier(&domain[..half].to_vec(), primitive_root, root_order);
    let right_zerofier = fast_zerofier(&domain[half..].to_vec(), primitive_root, root_order);

    let left_offset = fast_evaluate(
        &right_zerofier,
        &domain[..half].to_vec(),
        primitive_root,
        root_order,
    );
    let right_offset = fast_evaluate(
        &left_zerofier,
        &domain[half..].to_vec(),
        primitive_root,
        root_order,
    );

    let left_targets: Vec<_> = values[..half]
        .iter()
        .zip(left_offset.iter())
        .map(|(n, d)| n.div_ref(d))
        .collect();
    let right_targets: Vec<_> = values[half..]
        .iter()
        .zip(right_offset.iter())
        .map(|(n, d)| n.div_ref(d))
        .collect();

    let left_interpolant = fast_interpolate(
        &domain[..half].to_vec(),
        &left_targets,
        primitive_root,
        root_order,
    );
    let right_interpolant = fast_interpolate(
        &domain[half..].to_vec(),
        &right_targets,
        primitive_root,
        root_order,
    );

    left_interpolant.reduce() * right_zerofier.reduce()
        + right_interpolant.reduce() * left_zerofier.reduce()
}

pub fn fast_coset_evaluate<F>(
    polynomial: &Polynomial<F>,
    offset: &F,
    generator: &F,
    order: usize,
) -> Vec<F>
where
    F: Field,
{
    let scaled_polynomial = polynomial.scale(offset);
    let mut coefs = scaled_polynomial.coef.clone();
    for _ in 0..(order - polynomial.coef.len()) {
        coefs.push(F::zero());
    }
    ntt(&generator, &coefs)
}

pub fn fast_coset_divide<F>(
    lhs: &Polynomial<F>,
    rhs: &Polynomial<F>,
    offset: &F,
    primitive_root: &F,
    root_order: usize,
) -> Polynomial<F>
where
    F: Field,
{
    assert!(primitive_root.pow(root_order).is_one());
    assert!(!primitive_root.pow(root_order / 2).is_one());
    assert!(!rhs.is_zero());
    assert!(rhs.degree() < lhs.degree());

    if lhs.is_zero() {
        return Polynomial::zero();
    }

    let mut root = primitive_root.clone();
    let mut order = root_order.clone();
    let degree = std::cmp::max(lhs.degree(), rhs.degree());

    if degree < 8 {
        return lhs / rhs;
    }

    while degree < (order / 2).try_into().unwrap() {
        root = root.mul_ref(&root);
        order = order / 2;
    }

    let scaled_lhs = lhs.scale(offset);
    let scaled_rhs = rhs.scale(offset);

    let mut lhs_coefficients: Vec<_> = scaled_lhs.coef[..lhs.degree() as usize + 1].to_vec();
    while lhs_coefficients.len() < order {
        lhs_coefficients.push(F::zero());
    }
    let mut rhs_coefficients: Vec<_> = scaled_rhs.coef[..rhs.degree() as usize + 1].to_vec();
    while rhs_coefficients.len() < order {
        rhs_coefficients.push(F::zero());
    }

    let lhs_codeword = ntt(&root, &lhs_coefficients);
    let rhs_codeword = ntt(&root, &rhs_coefficients);

    let quotient_codeword: Vec<_> = lhs_codeword
        .into_iter()
        .zip(rhs_codeword.iter())
        .map(|(el, r)| el.div_ref(r))
        .collect();
    let scaled_quotient_coefficients = intt(&root, &quotient_codeword);
    let scaled_quotient = Polynomial {
        coef: scaled_quotient_coefficients[..(lhs.degree() as usize - rhs.degree() as usize + 1)]
            .to_vec(),
    };

    scaled_quotient.scale(&offset.inverse())
}

#[cfg(test)]
mod tests {
    use std::primitive;

    use num_bigint::BigInt;

    use super::*;

    use crate::define_extension_field;
    use crate::modules::algebra::field::FiniteFieldElement;
    use crate::modules::algebra::ring::Ring;
    use crate::modules::zkstark::fri::{get_nth_root_of_m128, M128};
    // 1 + 407 * (1 << 119)

    #[test]
    fn test_ntt() {
        let logn = 8;
        let n = 1 << logn;
        let primitive_root = get_nth_root_of_m128(&BigInt::from(n));

        let coef: Vec<FiniteFieldElement<M128>> = (0..n)
            .into_iter()
            .map(|i| FiniteFieldElement::<M128>::from_value(i + 1))
            .collect();
        let poly = Polynomial { coef: coef.clone() };

        let values = ntt(&primitive_root, &coef);
        let values_again = poly.eval_domain(
            &(0..values.len())
                .map(|i| primitive_root.pow(i))
                .collect::<Vec<_>>(),
        );
        assert_eq!(values.len(), values_again.len());
        for i in 0..values.len() {
            assert_eq!(values[i], values_again[i]);
        }

        let coef_again = intt(&primitive_root, &values);
        assert_eq!(coef.len(), coef_again.len());
        for i in 0..coef.len() {
            assert_eq!(coef[i], coef_again[i]);
        }
    }
}
