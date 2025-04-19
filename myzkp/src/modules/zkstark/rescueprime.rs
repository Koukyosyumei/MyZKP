use std::collections::HashMap;
use std::hash::Hash;
use std::thread::current;

use blake2::{digest::consts::U32, Blake2b, Digest};
use num_traits::{zero, One, Zero};

use crate::modules::algebra::field::{Field, FiniteFieldElement, ModulusValue};
use crate::modules::algebra::merkle::Merkle;
use crate::modules::algebra::mpolynomials::MPolynomial;
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::ring::Ring;
use crate::modules::zkstark::fiat_shamir::FiatShamirTransformer;
use crate::modules::zkstark::fri::FRI;
use crate::modules::zkstark::stark::Stark;

pub struct RescuePrime<M: ModulusValue> {
    pub p: usize,
    pub field: FiniteFieldElement<M>,
    pub m: usize,
    pub rate: usize,
    pub capacity: usize,
    pub n: usize,
    pub alpha: usize,
    pub alphainv: usize,
    pub mds: Vec<Vec<FiniteFieldElement<M>>>,
    pub mdsinv: Vec<Vec<FiniteFieldElement<M>>>,
    pub round_constants: Vec<FiniteFieldElement<M>>,
}

impl<M: ModulusValue> RescuePrime<M> {
    pub fn hash(&self, input_element: &FiniteFieldElement<M>) {
        // absorb
        let mut state = vec![input_element.clone()];
        for _ in 0..(self.m - 1) {
            state.push(FiniteFieldElement::<M>::zero());
        }

        // permutation
        for r in 0..self.n {
            // forward half-round
            // s-box
            for i in 0..self.m {
                state[i] = state[i].pow(self.alpha);
            }
            // matrix
            let mut temp: Vec<_> = (0..self.m)
                .into_iter()
                .map(|_| FiniteFieldElement::<M>::zero())
                .collect();
            for i in 0..self.m {
                for j in 0..self.m {
                    temp[i] = temp[i].clone() + self.mds[i][j].clone() * state[j].clone();
                }
            }
            // constants
            for i in 0..self.m {
                state[i] = temp[i].clone() + self.round_constants[2 * r * self.m + i].clone();
            }
        }
    }
}
