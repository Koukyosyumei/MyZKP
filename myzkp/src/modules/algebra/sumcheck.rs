use std::collections::HashMap;

use num_traits::{One, Zero};

use crate::modules::algebra::curve::bn128::FqOrder;
use crate::modules::algebra::fiat_shamir::FiatShamirTransformer;
use crate::modules::algebra::field::Field;
use crate::modules::algebra::gemini::{
    commit_gemini, open_gemini, split_and_fold, verify_gemini, CommitmentGemini, ProofGemini,
};
use crate::modules::algebra::kzg::PublicKeyKZG;
use crate::modules::algebra::mpolynomials::MPolynomial;
use crate::modules::algebra::polynomial::Polynomial;
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
            let bit = (n >> i) & 1;
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
                .unwrap_or(&F::zero())
                .clone()
        })
        .collect::<Vec<_>>()
}

pub struct SumCheckProof {
    pub h: FqOrder,
    pub el: usize,
    pub gs: Vec<MPolynomial<FqOrder>>,
    pub c_g: CommitmentGemini,
    pub pi: ProofGemini,
}

pub fn commit_sumcheck(
    g: &MPolynomial<FqOrder>,
    rs: &Vec<FqOrder>,
    pk: &PublicKeyKZG,
) -> (CommitmentGemini, Vec<Polynomial<FqOrder>>) {
    let coefs = get_coefs_in_order(&g);
    let fs = split_and_fold(&coefs, &rs).unwrap();
    (commit_gemini(&fs, &pk), fs)
}

pub fn prove_sumcheck(g: &MPolynomial<FqOrder>, h: &FqOrder, pk: &PublicKeyKZG) -> SumCheckProof {
    let mut proof_stream = FiatShamirTransformer::new();

    let el = g.get_num_vars();
    proof_stream.push(&vec![bincode::serialize(&el).expect("Serialization failed")]);
    proof_stream.push(&vec![bincode::serialize(&h).expect("Serialization failed")]);

    let mut gs = vec![];
    let mut rs = vec![];

    let g_0 = build_gj_from_prefix(&g, &vec![]);
    proof_stream.push(&vec![
        bincode::serialize(&g_0).expect("Serialization failed")
    ]);
    gs.push(g_0);
    let r_0 = FqOrder::sample(&proof_stream.prover_fiat_shamir(32));
    rs.push(r_0);

    for _ in 1..el {
        let g_j = build_gj_from_prefix(&g, &rs);
        proof_stream.push(&vec![
            bincode::serialize(&g_j).expect("Serialization failed")
        ]);
        let r_j = FqOrder::sample(&proof_stream.prover_fiat_shamir(32));
        rs.push(r_j);
        gs.push(g_j);
    }

    let beta = FqOrder::sample(&proof_stream.prover_fiat_shamir(32));
    let (c_g, fs) = commit_sumcheck(g, &rs, pk);
    let pi = open_gemini(&fs, &beta, &pk);

    SumCheckProof {
        h: h.clone(),
        el,
        gs,
        c_g,
        pi,
    }
}

pub fn verify_sumcheck(proof: &SumCheckProof, pk: &PublicKeyKZG) -> bool {
    let mut proof_stream = FiatShamirTransformer::new();
    proof_stream.push(&vec![
        bincode::serialize(&proof.el).expect("Serialization failed")
    ]);
    proof_stream.push(&vec![
        bincode::serialize(&proof.h).expect("Serialization failed")
    ]);

    if proof.h != sumcheck_fold(&proof.gs[0], 0) {
        return false;
    }

    let mut rs = vec![];
    proof_stream.push(&vec![
        bincode::serialize(&proof.gs[0]).expect("Serialization failed")
    ]);
    let r_0 = FqOrder::sample(&proof_stream.prover_fiat_shamir(32));
    rs.push(r_0);

    for j in 1..proof.el {
        let mut rs_j_minus_1_vec = vec![FqOrder::zero(); proof.el];
        rs_j_minus_1_vec[j - 1] = rs[j - 1].clone();
        if proof.gs[j - 1].evaluate(&rs_j_minus_1_vec) != sumcheck_fold(&proof.gs[j], j) {
            return false;
        }
        proof_stream.push(&vec![
            bincode::serialize(&proof.gs[j]).expect("Serialization failed")
        ]);
        let r_j = FqOrder::sample(&proof_stream.prover_fiat_shamir(32));
        rs.push(r_j);
    }

    let beta = FqOrder::sample(&proof_stream.prover_fiat_shamir(32));
    let mut rs_el_minus_1_vec = vec![FqOrder::zero(); proof.el];
    rs_el_minus_1_vec[proof.el - 1] = rs[proof.el - 1].clone();
    let g_el_minus_1_at_el_minus_1 = proof.gs[proof.el - 1].evaluate(&rs_el_minus_1_vec);

    verify_gemini(
        &rs,
        &g_el_minus_1_at_el_minus_1,
        &beta,
        &proof.c_g,
        &proof.pi,
        &pk,
    )
}

#[cfg(test)]
mod tests {
    use std::hash::Hash;

    use num_traits::One;

    use super::*;

    use crate::modules::algebra::curve::bn128::BN128;
    use crate::modules::algebra::field::{FiniteFieldElement, ModEIP197};
    use crate::modules::algebra::kzg::setup_kzg_with_full_g2;
    use crate::modules::algebra::ring::Ring;

    #[test]
    fn test_sumcheck_pipeline() {
        let g1 = BN128::generator_g1();
        let g2 = BN128::generator_g2();
        let pk = setup_kzg_with_full_g2(&g1, &g2, 8);

        type F = FiniteFieldElement<ModEIP197>;
        let mut dict = HashMap::new();
        dict.insert(vec![0, 0, 0], F::from_value(1));
        dict.insert(vec![1, 0, 0], F::from_value(2));
        dict.insert(vec![0, 1, 0], F::from_value(3));
        dict.insert(vec![0, 1, 1], F::from_value(4));
        dict.insert(vec![1, 1, 1], F::from_value(5));
        let g = MPolynomial::new(dict);
        let h = sum_over_boolean_hypercube(&g);

        let proof = prove_sumcheck(&g, &h, &pk);
        assert!(verify_sumcheck(&proof, &pk));
    }

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
    fn test_first_round() {
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
