use std::fmt;

use num_traits::Zero;

use crate::modules::algebra::curve::bn128::{optimal_ate_pairing, FqOrder, G1Point, G2Point};
use crate::modules::algebra::field::Field;
use crate::modules::algebra::kzg::{
    batch_open_kzg, batch_verify_kzg, commit_kzg, open_kzg, prove_degree_bound, setup_kzg,
    verify_degree_bound, verify_kzg, BatchProofKZG, CommitmentKZG, ProofDegreeBound, ProofKZG,
    PublicKeyKZG,
};
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
            let bit = (n >> (self.len - 1 - i)) & 1;
            combination.push(bit as u8);
        }

        self.current += 1;

        Some(combination)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bitcombinations() {
        let combinations = BitCombinations::new(3, 0);
        let vs = combinations.into_iter().collect::<Vec<_>>();
        assert_eq!(vs.len(), 8);

        let combinations = BitCombinations::new(3, 1);
        let vs = combinations.into_iter().collect::<Vec<_>>();
        assert_eq!(vs.len(), 7);
    }
}
