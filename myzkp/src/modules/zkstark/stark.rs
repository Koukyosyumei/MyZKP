use blake2::{digest::consts::U32, Blake2b, Digest};
use num_traits::One;

use crate::modules::algebra::field::{Field, FiniteFieldElement, ModulusValue};
use crate::modules::algebra::merkle::Merkle;
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::ring::Ring;
use crate::modules::zkstark::fiat_shamir::FiatShamirTransformer;
use crate::modules::zkstark::fri::FRI;

pub struct Stark<M: ModulusValue> {
    pub expansion_factor: usize,
    pub num_colinearity_checks: usize,
    pub security_level: usize,
    pub num_randomizers: usize,
    pub num_registers: usize,
    pub original_trace_length: usize,
    pub omicron: FiniteFieldElement<M>,
    pub omicron_domain: Vec<FiniteFieldElement<M>>,
    pub fri: FRI<M>,
}

pub type Trace<M> = Vec<Vec<FiniteFieldElement<M>>>;
pub type Boundary<M> = Vec<(usize, usize, FiniteFieldElement<M>)>;

impl<M: ModulusValue> Stark<M> {
    fn transition_zerofier(&self) -> Polynomial<FiniteFieldElement<M>> {
        let domain = &self.omicron_domain[0..(self.original_trace_length - 1)];
        Polynomial::from_monomials(domain)
    }

    fn boundary_zerofiers(&self, boundary: &Boundary<M>) -> Vec<Polynomial<FiniteFieldElement<M>>> {
        let mut zerofiers = Vec::new();
        for s in 0..self.num_registers {
            let mut points = Vec::new();
            for (c, r, v) in boundary {
                if *r == s {
                    points.push(self.omicron.pow(*c));
                }
            }
            zerofiers.push(Polynomial::from_monomials(&points));
        }
        zerofiers
    }

    fn boundary_interpolants(
        &self,
        boundary: &Boundary<M>,
    ) -> Vec<Polynomial<FiniteFieldElement<M>>> {
        let mut interpolants = Vec::new();
        for s in 0..self.num_registers {
            let mut points = Vec::new();
            let mut domain = Vec::new();
            let mut values = Vec::new();
            for (c, r, v) in boundary {
                if *r == s {
                    points.push((c, v));
                    domain.push(self.omicron.pow(*c));
                    values.push(v.clone());
                }
            }
            interpolants.push(Polynomial::interpolate(&domain, &values));
        }
        interpolants
    }

    fn boundary_quotient_degree_bounds(
        &self,
        randomized_trace_length: usize,
        boundary: &Boundary<M>,
    ) -> Vec<usize> {
        let randomized_trace_degree = randomized_trace_length - 1;
        self.boundary_zerofiers(boundary)
            .iter()
            .map(|bz| randomized_trace_degree - bz.degree() as usize)
            .collect()
    }

    pub fn prove(
        &self,
        trace: &mut Trace<M>,
        boundary: &Boundary<M>,
        proof_stream: &mut FiatShamirTransformer,
    ) {
        // concatenate randomizers
        for _ in 0..self.num_randomizers {
            trace.push(
                (0..self.num_registers)
                    .into_iter()
                    .map(|_| FiniteFieldElement::<M>::random_element(&[]))
                    .collect(),
            );
        }

        // interpolate
        let trace_domain: Vec<_> = (0..trace.len())
            .into_iter()
            .map(|i| self.omicron.pow(i))
            .collect();
        let mut trace_polynomials = Vec::new();
        for s in 0..self.num_registers {
            let single_trace: Vec<_> = (0..trace.len())
                .into_iter()
                .map(|c| trace[c][s].clone())
                .collect();
            trace_polynomials.push(Polynomial::interpolate(&trace_domain, &single_trace));
        }

        // subtract boundary interpolants and divide out boundary zerofiers
        let mut boundary_quotients = Vec::new();
        for s in 0..self.num_registers {
            let interpolant = self.boundary_interpolants(&boundary)[s].clone();
            let zerofier = self.boundary_zerofiers(&boundary)[s].clone();
            let quotient = (trace_polynomials[s].clone() - interpolant) / zerofier;
            boundary_quotients.push(quotient);
        }

        // commit to boundary quotients
        let fri_domain = self.fri.eval_domain();
        let mut boundary_quotient_codewords = Vec::new();
        //let mut boundary_quotient_Merkle_roots = Vec::new();
        for s in 0..self.num_registers {
            boundary_quotient_codewords.push(boundary_quotients[s].eval_domain(&fri_domain));
            let tmp: Vec<_> = boundary_quotient_codewords[s]
                .iter()
                .map(|c| bincode::serialize(c).expect("Serialization failed"))
                .collect();
            let merkle_root = Merkle::commit(&tmp);
            proof_stream.push(&vec![merkle_root]);
        }

        // symbolically evaluate transition constraints
    }
}
