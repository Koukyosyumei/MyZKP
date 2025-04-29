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
use crate::modules::zkstark::stark::{Boundary, Stark, TransitionConstraints};

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
    pub fn hash(&self, input_element: &FiniteFieldElement<M>) -> FiniteFieldElement<M> {
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

            // backward half-round
            // s-box
            for i in 0..self.m {
                state[i] = state[i].pow(self.alphainv);
            }
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
                state[i] =
                    temp[i].clone() + self.round_constants[2 * r * self.m + self.m + i].clone();
            }
        }

        state[0].clone()
    }

    fn round_constants_polynomials(
        &self,
        omicron: &FiniteFieldElement<M>,
    ) -> (
        Vec<MPolynomial<FiniteFieldElement<M>>>,
        Vec<MPolynomial<FiniteFieldElement<M>>>,
    ) {
        let mut first_step_constants = Vec::new();
        for i in 0..self.m {
            let domain: Vec<_> = (0..self.n).map(|r| omicron.pow(r)).collect();
            let values: Vec<_> = (0..self.n)
                .map(|r| self.round_constants[2 * r * self.m + i].clone())
                .collect();
            let univariate = Polynomial::interpolate(&domain, &values);
            let multivariate = MPolynomial::lift(&univariate, 0);
            first_step_constants.push(multivariate);
        }

        let mut second_step_constants = Vec::new();
        for i in 0..self.m {
            let domain: Vec<_> = (0..self.n).map(|r| omicron.pow(r)).collect();
            let values: Vec<_> = (0..self.n)
                .map(|r| self.round_constants[2 * r * self.m + self.m + i].clone())
                .collect();
            let univariate = Polynomial::interpolate(&domain, &values);
            let multivariate = MPolynomial::lift(&univariate, 0);
            second_step_constants.push(multivariate);
        }

        (first_step_constants, second_step_constants)
    }

    pub fn transition_constraints(
        &self,
        omicron: FiniteFieldElement<M>,
    ) -> TransitionConstraints<M> {
        let (first_step_constants, second_step_constants) =
            self.round_constants_polynomials(&omicron);

        let variables = MPolynomial::<FiniteFieldElement<M>>::variables(1 + 2 * self.m);
        //let cycle_index = &variables[0];
        let previous_state = &variables[1..(1 + self.m)];
        let next_state = &variables[(1 + self.m)..(1 + 2 * self.m)];

        let mut air = Vec::new();
        for i in 0..self.m {
            let mut lhs = MPolynomial::constant(FiniteFieldElement::<M>::zero());
            for k in 0..self.m {
                lhs = lhs
                    + MPolynomial::constant(self.mds[i][k].clone())
                        * (previous_state[k].pow(self.alpha));
            }
            lhs = lhs + first_step_constants[i].clone();

            let mut rhs = MPolynomial::constant(FiniteFieldElement::<M>::zero());
            for k in 0..self.m {
                rhs = rhs
                    + MPolynomial::constant(self.mdsinv[i][k].clone())
                        * (&next_state[k] - &second_step_constants[k]);
            }
            rhs = rhs.pow(self.alpha);
            air.push(lhs - rhs);
        }

        air
    }

    pub fn boundary_constraints(&self, output_element: FiniteFieldElement<M>) -> Boundary<M> {
        let mut constraints = Vec::new();
        constraints.push((0, 1, FiniteFieldElement::<M>::zero()));
        constraints.push((self.n, 0, output_element));
        constraints
    }
}
