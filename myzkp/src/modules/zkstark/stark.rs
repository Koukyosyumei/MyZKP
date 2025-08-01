use std::collections::HashMap;
use std::str::FromStr;

use blake2::{digest::consts::U32, Blake2b, Digest};
use num_bigint::BigInt;
use num_traits::{One, Zero};

use crate::modules::algebra::field::{Field, FiniteFieldElement, ModulusValue};
use crate::modules::algebra::merkle::{Merkle, MerklePath, MerkleRoot};
use crate::modules::algebra::mpolynomials::MPolynomial;
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::ring::Ring;
use crate::modules::zkstark::fiat_shamir::FiatShamirTransformer;
use crate::modules::zkstark::fri::{FriProof, FRI};

use super::fri::{get_nth_root_of_m128, M128};

pub type Trace<F> = Vec<Vec<F>>;
pub type Boundary<F> = Vec<(usize, usize, F)>;
pub type TransitionConstraints<F> = Vec<MPolynomial<F>>;
pub struct StarkProof<F: Field> {
    pub fri_proof: FriProof<F>,
    pub bqc_roots: Vec<MerkleRoot>,
    pub bqc_points: Vec<F>,
    pub bqc_paths: Vec<MerklePath>,
    pub rdc_root: MerkleRoot,
    pub rdc_points: Vec<F>,
    pub rdc_paths: Vec<MerklePath>,
}

pub struct Stark<F: Field> {
    pub expansion_factor: usize,
    pub num_colinearity_checks: usize,
    pub security_level: usize,
    pub num_randomizers: usize,
    pub num_registers: usize,
    pub original_trace_length: usize,
    pub generator: F,
    pub omega: F,
    pub omicron: F,
    pub omicron_domain: Vec<F>,
    pub fri: FRI<F>,
}

impl<F: Field> Stark<F> {
    fn transition_degree_bounds(
        &self,
        transition_constraints: &TransitionConstraints<F>,
    ) -> Vec<usize> {
        let mut point_degrees = vec![1];
        for _ in 0..(2 * self.num_registers) {
            point_degrees.push(self.original_trace_length + self.num_randomizers - 1);
        }
        transition_constraints
            .iter()
            .map(|a| {
                a.dictionary
                    .iter()
                    .map(|(k, v)| point_degrees.iter().zip(k).map(|(r, el)| r * el).sum())
                    .max()
                    .unwrap_or(0)
            })
            .collect()
    }

    fn transition_quotient_degree_bounds(
        &self,
        transition_constraints: &Vec<MPolynomial<F>>,
    ) -> Vec<usize> {
        self.transition_degree_bounds(transition_constraints)
            .iter()
            .map(|d| d - (self.original_trace_length - 1))
            .collect()
    }

    fn max_degree(&self, transition_constraints: &Vec<MPolynomial<F>>) -> usize {
        let binding = self.transition_quotient_degree_bounds(transition_constraints);
        let md = binding.iter().max().unwrap();
        (1 << (format!("{:b}", md).len())) - 1
    }

    fn transition_zerofier(&self) -> Polynomial<F> {
        let domain = &self.omicron_domain[0..(self.original_trace_length - 1)];
        Polynomial::from_monomials(domain)
    }

    fn boundary_zerofiers(&self, boundary: &Boundary<F>) -> Vec<Polynomial<F>> {
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

    fn boundary_interpolants(&self, boundary: &Boundary<F>) -> Vec<Polynomial<F>> {
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
        boundary: &Boundary<F>,
    ) -> Vec<usize> {
        let randomized_trace_degree = randomized_trace_length - 1;
        self.boundary_zerofiers(boundary)
            .iter()
            .map(|bz| randomized_trace_degree - bz.degree() as usize)
            .collect()
    }

    fn sample_weights(&self, number: usize, randomness: Vec<u8>) -> Vec<F> {
        let mut result = Vec::new();
        for i in 0..number {
            let mut tmp_r = randomness.clone();
            let mut tmp_i: Vec<u8> = i.to_le_bytes().to_vec();
            tmp_r.append(&mut tmp_i);

            let mut hasher = Blake2b::<U32>::new();
            hasher.update(tmp_r);
            let hash_result = hasher.finalize();
            result.push(F::sample(&hash_result));
        }
        result
    }

    pub fn prove(
        &self,
        trace: &mut Trace<F>,
        boundary: &Boundary<F>,
        transition_constraints: &TransitionConstraints<F>,
    ) -> StarkProof<F> {
        let mut proof_stream = FiatShamirTransformer::new();

        // concatenate randomizers
        for _ in 0..self.num_randomizers {
            trace.push(
                (0..self.num_registers)
                    .into_iter()
                    .map(|_| F::random_element(&[]))
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
        let mut bqc_roots = Vec::new();
        for s in 0..self.num_registers {
            boundary_quotient_codewords.push(boundary_quotients[s].eval_domain(&fri_domain));
            let tmp: Vec<_> = boundary_quotient_codewords[s]
                .iter()
                .map(|c| bincode::serialize(c).expect("Serialization failed"))
                .collect();
            let merkle_root = Merkle::commit(&tmp);
            bqc_roots.push(merkle_root.clone());
            proof_stream.push(&vec![merkle_root]);
        }

        // symbolically evaluate transition constraints
        let mut point = vec![Polynomial {
            coef: vec![F::zero(), F::one()],
        }];
        for tp in &trace_polynomials {
            point.push(tp.clone());
        }
        for tp in trace_polynomials {
            point.push(tp.scale(&self.omicron));
        }
        let transition_polynomials: Vec<_> = transition_constraints
            .iter()
            .map(|a| a.evaluate_symbolic(&point))
            .collect();

        // divide out zerofier
        let transition_quotients: Vec<_> = transition_polynomials
            .iter()
            .map(|tp| tp / &self.transition_zerofier())
            .collect();

        // commit to randomizer polynomial
        let randomizer_polynomial = Polynomial {
            coef: (0..(self.max_degree(transition_constraints) + 1))
                .into_iter()
                .map(|_| F::random_element(&[]))
                .collect(),
        };
        let randomizer_codeword = randomizer_polynomial.eval_domain(&fri_domain);
        let tmp: Vec<_> = randomizer_codeword
            .iter()
            .map(|c| bincode::serialize(c).expect("Serialization failed"))
            .collect();
        let randomizer_root = Merkle::commit(&tmp);
        proof_stream.push(&vec![randomizer_root.clone()]);

        // get weights for nonlinear combination
        let weights = self.sample_weights(
            1 + 2 * transition_quotients.len() + 2 * boundary_quotients.len(),
            proof_stream.prover_fiat_shamir(32),
        );

        // compute terms of nonlinear combination polynomial
        let x = Polynomial {
            coef: vec![F::zero(), F::one()],
        };
        let mut terms = vec![randomizer_polynomial];
        for i in 0..transition_quotients.len() {
            terms.push(transition_quotients[i].clone());
            let shift = self.max_degree(transition_constraints)
                - self.transition_quotient_degree_bounds(transition_constraints)[i];
            terms.push(x.pow(shift) * transition_quotients[i].clone());
        }
        for i in 0..self.num_registers {
            terms.push(boundary_quotients[i].clone());
            let shift = self.max_degree(transition_constraints)
                - self.boundary_quotient_degree_bounds(trace.len(), boundary)[i];
            terms.push(x.pow(shift) * boundary_quotients[i].clone());
        }

        // take weighted sum
        let mut combination = Polynomial::<F>::zero();
        for i in 0..terms.len() {
            let tmp = Polynomial {
                coef: vec![weights[i].clone()],
            };
            combination = combination + tmp * terms[i].clone();
        }

        // compute matching codeword
        let combined_codeword = combination.eval_domain(&fri_domain);

        // prove low degree of combination polynomial
        let mut fri_proof = self.fri.prove(&combined_codeword);
        fri_proof.top_level_indices.sort();
        let mut duplicated_indices = fri_proof.top_level_indices.clone();
        for i in &fri_proof.top_level_indices {
            duplicated_indices.push((i + self.expansion_factor) % self.fri.domain_length);
        }
        for i in duplicated_indices.clone() {
            duplicated_indices.push((i + (self.fri.domain_length / 2)) % self.fri.domain_length);
        }
        duplicated_indices.sort();

        // open indicated positions in the boundary quotient codewords
        let mut bqc_points = Vec::new();
        let mut bqc_paths = Vec::new();
        for bqc in boundary_quotient_codewords {
            let serialized_bqc: Vec<_> = bqc
                .iter()
                .map(|b| bincode::serialize(&b).expect("Serialization failed"))
                .collect();
            for i in &duplicated_indices {
                bqc_points.push(bqc[*i].clone());
                bqc_paths.push(Merkle::open(*i, &serialized_bqc));
            }
        }

        // as well as in the randomizer
        let mut rdc_points = Vec::new();
        let mut rdc_paths = Vec::new();
        let serialized_randomizer_codeword: Vec<_> = randomizer_codeword
            .iter()
            .map(|b| bincode::serialize(&b).expect("Serialization failed"))
            .collect();
        for i in &duplicated_indices {
            rdc_points.push(randomizer_codeword[*i].clone());
            rdc_paths.push(Merkle::open(*i, &serialized_randomizer_codeword));
        }

        StarkProof {
            fri_proof: fri_proof,
            bqc_roots: bqc_roots,
            bqc_points: bqc_points,
            bqc_paths: bqc_paths,
            rdc_root: randomizer_root,
            rdc_points: rdc_points,
            rdc_paths: rdc_paths,
        }
    }

    pub fn verify(
        &self,
        proof: &StarkProof<F>,
        transition_constraints: &TransitionConstraints<F>,
        boundary: &Boundary<F>,
    ) -> bool {
        let mut proof_stream = FiatShamirTransformer::new();

        // infer trace length from boundary conditions
        let original_trace_length = 1 + boundary.iter().map(|(c, _, _)| c).max().unwrap();
        let randomized_trace_length = original_trace_length + self.num_randomizers;

        // get Merkle roots of boundary quotient codewords
        let boundary_quotient_roots = &proof.bqc_roots;
        for bqr in boundary_quotient_roots {
            proof_stream.push(&vec![bqr.clone()]);
        }

        // get Merkle root of randomizer polynomial
        let randomizer_root = &proof.rdc_root;
        proof_stream.push(&vec![randomizer_root.clone()]);

        // get weights for nonlinear combination
        let weights = self.sample_weights(
            1 + 2 * transition_constraints.len() + 2 * self.boundary_interpolants(boundary).len(),
            proof_stream.prover_fiat_shamir(32),
        );

        // verify low degree of combination polynomial
        let mut polynomial_values = Vec::new();
        let mut verifier_accepts = self.fri.verify(&proof.fri_proof, &mut polynomial_values);
        polynomial_values.sort_by_key(|iv| iv.0);
        if !verifier_accepts {
            return false;
        }

        let indices: Vec<_> = polynomial_values.iter().map(|(i, _)| i.clone()).collect();
        let values: Vec<_> = polynomial_values.iter().map(|(_, v)| v.clone()).collect();

        // read and verify leafs, which are elements of boundary quotient codewords
        let mut duplicated_indices = indices.clone();
        for i in &indices {
            duplicated_indices.push((i + self.expansion_factor) % self.fri.domain_length);
        }
        duplicated_indices.sort();

        let mut leafs = Vec::new();
        let mut ctr = 0;
        for r in 0..boundary_quotient_roots.len() {
            let mut tmp: HashMap<&usize, Vec<u8>> = HashMap::new();
            for i in &duplicated_indices {
                tmp.insert(
                    i,
                    bincode::serialize(&proof.bqc_points[ctr]).expect("Serialization failed"),
                );
                let path = &proof.bqc_paths[ctr];
                let verifier_accepts = verifier_accepts
                    && Merkle::verify(&boundary_quotient_roots[r], *i, &path, &tmp[&i]);
                if !verifier_accepts {
                    return false;
                }
                ctr += 1;
            }
            leafs.push(tmp);
        }

        // read and verify randomizer leafs
        let mut randomizer = HashMap::new();
        for (ctr, i) in duplicated_indices.iter().enumerate() {
            randomizer.insert(
                i,
                bincode::serialize(&proof.rdc_points[ctr]).expect("Serialization failed"),
            );
            let path = &proof.rdc_paths[ctr];
            verifier_accepts =
                verifier_accepts && Merkle::verify(&randomizer_root, *i, &path, &randomizer[&i]);
        }

        // verify leafs of combination polynomial
        for i in 0..indices.len() {
            let current_index = indices[i];

            // get trace values by applying a correction to the boundary quotient values (which are the leafs)
            let domain_current_index = self.generator.clone() * (self.omega.pow(current_index));
            let next_index = (current_index + self.expansion_factor) % self.fri.domain_length;
            let domain_next_index = self.generator.clone() * (self.omega.pow(next_index));
            let mut current_trace: Vec<_> = (0..self.num_registers)
                .into_iter()
                .map(|_| F::zero())
                .collect();
            let mut next_trace: Vec<_> = (0..self.num_registers)
                .into_iter()
                .map(|_| F::zero())
                .collect();
            for s in 0..self.num_registers {
                let zerofier = self.boundary_zerofiers(boundary)[s].clone();
                let interpolant = self.boundary_interpolants(boundary)[s].clone();
                let tmp_cur: F = bincode::deserialize(&leafs[s][&current_index])
                    .expect("Deserialization failed");
                let tmp_next: F =
                    bincode::deserialize(&leafs[s][&next_index]).expect("Deserialization failed");
                current_trace[s] = tmp_cur * zerofier.eval(&domain_current_index)
                    + interpolant.eval(&domain_current_index);
                next_trace[s] = tmp_next * zerofier.eval(&domain_next_index)
                    + interpolant.eval(&domain_next_index);
            }

            let mut point = vec![domain_current_index.clone()];
            point.append(&mut current_trace);
            point.append(&mut next_trace);
            let transition_constraints_values: Vec<_> = (0..transition_constraints.len())
                .into_iter()
                .map(|s| transition_constraints[s].evaluate(&point))
                .collect();

            // compute nonlinear combination
            let mut terms =
                vec![bincode::deserialize(&randomizer[&current_index])
                    .expect("Deserialization failed")];
            for s in 0..(transition_constraints_values.len()) {
                let tcv = transition_constraints_values[s].clone();
                let quotient = tcv / self.transition_zerofier().eval(&domain_current_index);
                terms.push(quotient.clone());
                let shift = self.max_degree(transition_constraints)
                    - self.transition_quotient_degree_bounds(transition_constraints)[s];
                terms.push(quotient * (domain_current_index.pow(shift)));
            }
            for s in 0..self.num_registers {
                let tmp = &leafs[s][&current_index];
                let bqv: F = bincode::deserialize(&tmp).expect("Deserialization failed");
                terms.push(bqv.clone());
                let shift = self.max_degree(transition_constraints)
                    - self.boundary_quotient_degree_bounds(randomized_trace_length, boundary)[s];
                terms.push(bqv * (domain_current_index.pow(shift)));
            }
            let mut combination = F::zero();
            for j in 0..terms.len() {
                combination = combination + (terms[j].mul_ref(&weights[j]));
            }

            // verify against combination polynomial value
            verifier_accepts = verifier_accepts && (combination == values[i]);
            if !verifier_accepts {
                return false;
            }
        }

        verifier_accepts
    }
}

pub fn initialize_stark_m128(
    expansion_factor: usize,
    num_colinearity_checks: usize,
    security_level: usize,
    num_registers: usize,
    num_cycles: usize,
    transition_constraints_degree: usize,
) -> Stark<FiniteFieldElement<M128>> {
    let generator = FiniteFieldElement::<M128>::from_value(
        BigInt::from_str("85408008396924667383611388730472331217").unwrap(),
    );
    let num_randomizers = 4 * num_colinearity_checks;
    let randomized_trace_length = num_cycles + num_randomizers;
    let omicron_domain_length = 1
        << (usize::BITS as usize
            - (randomized_trace_length * transition_constraints_degree).leading_zeros() as usize);
    let fri_domain_length = omicron_domain_length * expansion_factor;
    let omega = get_nth_root_of_m128(&BigInt::from(fri_domain_length));
    let omicron = get_nth_root_of_m128(&BigInt::from(omicron_domain_length));
    let omicron_domain: Vec<_> = (0..omicron_domain_length).map(|i| omicron.pow(i)).collect();
    let fri = FRI {
        offset: generator.clone(),
        omega: omega.clone(),
        domain_length: fri_domain_length,
        expansion_factor: expansion_factor,
        num_colinearity_tests: num_colinearity_checks,
    };

    Stark::<FiniteFieldElement<M128>> {
        expansion_factor: expansion_factor,
        num_colinearity_checks: num_colinearity_checks,
        security_level: security_level,
        num_randomizers: num_randomizers,
        num_registers: num_registers,
        original_trace_length: num_cycles,
        generator: generator,
        omega: omega,
        omicron: omicron,
        omicron_domain: omicron_domain,
        fri: fri,
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use num_bigint::BigInt;

    use super::*;

    use crate::modules::{
        algebra::field::FiniteFieldElement,
        zkstark::{fri::M128, rescueprime::RescuePrime},
    };

    // 1 + 407 * (1 << 119)

    #[test]
    fn test_stark() {
        let expansion_factor = 4;
        let num_colinearity_checks = 2;
        let security_level = 2;

        let rp = RescuePrime::new();
        let mut output_element =
            FiniteFieldElement::<M128>::from_value(BigInt::from_str("123456789").unwrap());

        for _ in 0..10 {
            let input_element = output_element.clone();
            output_element = rp.hash(&input_element);
            let num_cycles = rp.n + 1;
            let state_width = rp.m;

            let stark = initialize_stark_m128(
                expansion_factor,
                num_colinearity_checks,
                security_level,
                state_width,
                num_cycles,
                2,
            );

            let mut trace = rp.trace(&input_element);
            let air = rp.transition_constraints(&stark.omicron);
            let boundary = rp.boundary_constraints(&output_element);
            let proof = stark.prove(&mut trace, &boundary, &air);
            let result = stark.verify(&proof, &air, &boundary);
            assert!(result);

            let false_output_element = output_element.clone() + FiniteFieldElement::<M128>::one();
            let false_boundary = rp.boundary_constraints(&false_output_element);
            let false_proof = stark.prove(&mut trace, &false_boundary, &air);
            let false_result = stark.verify(&false_proof, &air, &boundary);
            assert!(!false_result);
        }
    }
}
