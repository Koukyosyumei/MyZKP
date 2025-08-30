use std::collections::HashMap;
use std::fmt;

use num_traits::Zero;

use crate::modules::algebra::curve::bn128::{optimal_ate_pairing, FqOrder, G1Point, G2Point};
use crate::modules::algebra::field::{Field, FiniteFieldElement, ModEIP197};
use crate::modules::algebra::gemini::{commit_gemini, open_gemini, split_and_fold, verify_gemini};
use crate::modules::algebra::kzg::{
    batch_open_kzg, batch_verify_kzg, commit_kzg, open_kzg, prove_degree_bound, setup_kzg,
    verify_degree_bound, verify_kzg, BatchProofKZG, CommitmentKZG, ProofDegreeBound, ProofKZG,
    PublicKeyKZG,
};
use crate::modules::algebra::mpolynomials::MPolynomial;
use crate::modules::algebra::ring::Ring;

pub struct BitCombinations {
    len: usize,
    current: usize,
}

impl BitCombinations {
    pub fn new(length: usize, current: usize) -> Self {
        assert!(
            length <= usize::BITS as usize,
            "Length exceeds available bits"
        );
        BitCombinations {
            len: length,
            current: current,
        }
    }
}

impl Iterator for BitCombinations {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        let total_combinations = 1_usize.checked_shl(self.len as u32)?;
        if self.current >= total_combinations {
            return None;
        }

        let n = self.current;
        let mut combination = Vec::with_capacity(self.len);

        for i in 0..self.len {
            let bit = (n >> (self.len - 1 - i)) & 1;
            combination.push(bit as u8);
        }

        self.current += 1;

        Some(combination)
    }
}

pub fn sum_over_boolean_hypercube<F: Field>(g: &MPolynomial<F>) -> F {
    let el = g.get_num_vars();
    let comb = BitCombinations::new(el, 0);
    let mut h = F::zero();
    for c in comb {
        let c_casted = c.iter().map(|v| F::from_value(*v)).collect::<Vec<_>>();
        h = h + g.evaluate(&c_casted);
    }
    h
}

pub fn build_gj_from_prefix<F: Field>(g: &MPolynomial<F>, rs: &Vec<F>) -> MPolynomial<F> {
    let el = g.get_num_vars();
    let j_sub_1 = rs.len();
    assert!(el >= 1 && el > j_sub_1, "invalid sizes for sum-check round");
    let comb = BitCombinations::new(el - 1 - j_sub_1, 0);

    let mut g_j = MPolynomial::zero();
    for c in comb {
        let mut h = HashMap::new();
        for (i, v) in rs.iter().enumerate() {
            h.insert(i, v.clone());
        }
        for (i, v) in c.iter().enumerate() {
            h.insert(i + 1 + j_sub_1, F::from_value(v.clone()));
        }
        g_j = &g_j + &(g.partial_evaluate(&h));
    }

    g_j
}

pub fn sumcheck_fold<F: Field>(g_j: &MPolynomial<F>, j: usize) -> F {
    let el = g_j.get_num_vars();
    let mut one = vec![F::zero(); el];
    one[j] = F::one();
    let zero = vec![F::zero(); el];
    g_j.evaluate(&one) + g_j.evaluate(&zero)
}

pub fn get_coefs_in_order<F: Field>(g: &MPolynomial<F>) -> Vec<F> {
    let el = g.get_num_vars();
    let comb = BitCombinations::new(el, 0);
    comb.into_iter()
        .map(|v| {
            g.dictionary
                .get(&v.iter().map(|u| *u as usize).collect::<Vec<_>>())
                .unwrap()
                .clone()
        })
        .collect::<Vec<_>>()
}

pub fn sumcheck(
    h: &FqOrder,
    g: &MPolynomial<FqOrder>,
    rs: &Vec<FqOrder>,
    pk: &PublicKeyKZG,
) -> bool {
    let el = g.get_num_vars();

    let coefs = get_coefs_in_order(&g);
    let fs = split_and_fold(&coefs, &rs).unwrap();
    let commitment = commit_gemini(&fs, &pk);

    let g_0 = build_gj_from_prefix(&g, &vec![]);
    if h != &sumcheck_fold(&g_0, 0) {
        return false;
    }

    let mut g_j = g_0.clone();
    let mut g_j_minus_1 = g_0.clone();
    for j in 1..el {
        g_j = build_gj_from_prefix(&g, &rs);
        let mut rs_j_minus_1_vec = vec![FqOrder::zero(); el];
        rs_j_minus_1_vec[j - 1] = rs[j - 1].clone();
        if g_j_minus_1.evaluate(&rs_j_minus_1_vec) != sumcheck_fold(&g_j, j) {
            return false;
        }
        g_j_minus_1 = g_j.clone();
    }

    let mut rs_el_minus_1_vec = vec![FqOrder::zero(); el];
    rs_el_minus_1_vec[el - 1] = rs[el - 1].clone();
    let g_el_minus_1_at_el_minus_1 = g_j.evaluate(&rs_el_minus_1_vec);

    let beta = FqOrder::random_element(&vec![]);
    let proof = open_gemini(&fs, &beta, &pk);
    verify_gemini(
        &rs,
        &g_el_minus_1_at_el_minus_1,
        &beta,
        &commitment,
        &proof,
        &pk,
    )
}

#[cfg(test)]
mod tests {
    use std::hash::Hash;

    use num_traits::One;

    use super::*;

    use crate::modules::algebra::field::{FiniteFieldElement, ModEIP197};
    use crate::modules::algebra::ring::Ring;

    #[test]
    fn test_bitcombinations() {
        let combinations = BitCombinations::new(3, 0);
        let vs = combinations.into_iter().collect::<Vec<_>>();
        assert_eq!(vs.len(), 8);

        let combinations = BitCombinations::new(3, 1);
        let vs = combinations.into_iter().collect::<Vec<_>>();
        assert_eq!(vs.len(), 7);
    }

    #[test]
    fn test_check() {
        type F = FiniteFieldElement<ModEIP197>;
        let mut dict = HashMap::new();
        dict.insert(vec![0, 0, 0], F::from_value(1));
        dict.insert(vec![1, 0, 0], F::from_value(2));
        dict.insert(vec![0, 1, 0], F::from_value(3));
        dict.insert(vec![0, 1, 1], F::from_value(4));
        dict.insert(vec![1, 1, 1], F::from_value(5));
        let g = MPolynomial::new(dict);

        let h = sum_over_boolean_hypercube(&g);
        assert_eq!(h, F::from_value(41));

        let g_0 = build_gj_from_prefix(&g, &vec![]);
        assert_eq!(h, sumcheck_fold(&g_0, 0));
    }
}
