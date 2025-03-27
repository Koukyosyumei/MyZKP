use blake2::{digest::consts::U32, Blake2b, Digest};
use num_traits::One;

use crate::modules::algebra::field::{Field, FiniteFieldElement, ModulusValue};
use crate::modules::algebra::merkle::Merkle;
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::ring::Ring;
use crate::modules::zkstark::fiat_shamir::FiatShamirTransformer;
use crate::modules::zkstark::fri::FRI;

struct Stark<M: ModulusValue> {
    pub expansion_factor: usize,
    pub num_colinearity_checks: usize,
    pub security_level: usize,
    pub num_randomizers: usize,
    pub num_registers: usize,
    pub original_trace_length: usize,
    pub omicron: FiniteFieldElement<M>,
}

pub type Trace<M> = Vec<Vec<FiniteFieldElement<M>>>;
pub type Boundary<M> = Vec<(usize, usize, FiniteFieldElement<M>)>;

impl<M: ModulusValue> Stark<M> {
    pub fn prove(&self, trace: Trace<M>) {}
}
