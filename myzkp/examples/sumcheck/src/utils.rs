use std::collections::HashMap;
use std::sync::Arc;

use num_bigint::{BigInt, Sign};
use num_traits::Zero;

use myzkp::modules::algebra::field::{Field, FiniteFieldElement, ModEIP197};
use myzkp::modules::algebra::mpolynomials::MPolynomial;
use myzkp::modules::algebra::polynomial::Polynomial;
use myzkp::modules::algebra::ring::Ring;
use myzkp::modules::algebra::sumcheck::sum_over_boolean_hypercube;

#[derive(Debug)]
pub struct BitCombinationsDictOrder {
    len: usize,
    current: usize,
}

impl BitCombinationsDictOrder {
    pub fn new(length: usize) -> Self {
        assert!(
            length <= usize::BITS as usize,
            "Length exceeds available bits"
        );
        BitCombinationsDictOrder {
            len: length,
            current: 0,
        }
    }
}

impl Iterator for BitCombinationsDictOrder {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        let total_combinations = 1_usize.checked_shl(self.len as u32)?;
        if self.current >= total_combinations {
            return None;
        }

        let n = self.current;
        let mut combination = Vec::with_capacity(self.len);
        for i in (0..self.len).rev() {
            let bit = (n >> i) & 1;
            combination.push(bit as u8);
        }

        self.current += 1;
        Some(combination)
    }
}

pub type F = FiniteFieldElement<ModEIP197>;

pub fn field_to_bytes(x: &F) -> Vec<u8> {
    let (_, mut digits) = x.value.to_u64_digits();
    digits.resize(4, 0);
    digits.iter().flat_map(|d| d.to_le_bytes()).collect()
}

pub fn fields_to_bytes(xs: &Vec<F>) -> Vec<u8> {
    xs.iter().flat_map(|x| field_to_bytes(x)).collect()
}

pub fn fields_from_bytes(bytes: &[u8]) -> Vec<F> {
    const FIELD_BYTES: usize = 32;

    bytes
        .chunks(FIELD_BYTES)
        .map(|chunk| {
            let mut arr = [0u8; FIELD_BYTES];
            arr[..chunk.len()].copy_from_slice(chunk);
            F::from_value(BigInt::from_bytes_le(Sign::Plus, &arr))
        })
        .collect()
}

pub fn evals_over_boolean_hypercube(f: &MPolynomial<F>, result: &mut Vec<F>) {
    let el = f.get_num_vars();
    let comb = BitCombinationsDictOrder::new(el);
    for c in comb {
        let c_casted = c.iter().map(|v| F::from_value(*v)).collect::<Vec<_>>();
        result.push(f.evaluate(&c_casted));
    }
}

pub fn fold_factors_pointwise_cpu(evals: &mut [F], domain_size: usize, num_factors: usize) {
    for i in 0..domain_size {
        let mut product = evals[i].clone();
        for j in 1..num_factors {
            product = product * evals[i + j * domain_size].clone();
        }
        evals[i] = product;
    }
}

/// Corresponds to the `fold_into_half` CUDA kernel.
///
/// Performs one round of folding by binding the current variable to a challenge.
/// The update rule is `p'(x) = p(0,x) + r * (p(1,x) - p(0,x))`.
pub fn fold_into_half_cpu(
    evals: &mut [F],
    num_factors: usize,
    domain_size: usize, // Original domain size for offset calculation
    num_remaining_vars: usize,
    challenge: F,
) {
    let stride = 1 << (num_remaining_vars - 1);
    for factor_idx in 0..num_factors {
        let offset = factor_idx * domain_size;
        for i in 0..stride {
            let idx0 = offset + i;
            let idx1 = offset + i + stride;
            let val0 = &evals[idx0];
            let val1 = &evals[idx1];

            let diff = val1.sub_ref(val0);
            evals[idx0] = val0.add_ref(&(&challenge.mul_ref(&diff)));
        }
    }
}

/// Corresponds to the `eval_folded_poly` CUDA kernel.
///
/// Evaluates a partially folded polynomial at a specific point for the current variable.
/// The evaluation rule is `p(c, x) = p(0,x) + c * (p(1,x) - p(0,x))`.
pub fn eval_folded_poly_cpu(
    evals: &[F],
    num_factors: usize,
    domain_size: usize, // Original domain size
    num_vars: usize,
    eval_point: &F,
) -> Vec<F> {
    let stride = 1 << (num_vars - 1);
    let mut result = vec![F::zero(); num_factors * stride];

    for factor_idx in 0..num_factors {
        let poly_offset = factor_idx * domain_size;
        let result_offset = factor_idx * stride;

        for i in 0..stride {
            let idx0 = poly_offset + i;
            let idx1 = poly_offset + i + stride;
            let val0 = &evals[idx0];
            let val1 = &evals[idx1];

            let diff = val1.sub_ref(val0);
            result[result_offset + i] = val0.add_ref(&eval_point.mul_ref(&diff));
        }
    }
    result
}

pub fn sum_cpu(data: &[F], len: usize) -> F {
    let mut result = F::zero();
    for i in 0..len {
        result = result.add_ref(&data[i]);
    }
    result
}
