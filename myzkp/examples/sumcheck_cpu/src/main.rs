use std::collections::HashMap;
use std::fs;
use std::hash::Hash;
use std::sync::Arc;

use num_bigint::{BigInt, Sign};
use num_traits::identities::One;
use num_traits::Zero;

use myzkp::modules::algebra::fiat_shamir::FiatShamirTransformer;
use myzkp::modules::algebra::field::Field;
use myzkp::modules::algebra::field::{FiniteFieldElement, ModEIP197};
use myzkp::modules::algebra::mpolynomials::MPolynomial;
use myzkp::modules::algebra::polynomial::Polynomial;
use myzkp::modules::algebra::ring::Ring;
use myzkp::modules::algebra::sumcheck::sum_over_boolean_hypercube;

#[derive(Debug)]
struct BitCombinationsDictOrder {
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

type F = FiniteFieldElement<ModEIP197>;

fn field_to_bytes(x: &F) -> Vec<u8> {
    let (_, mut digits) = x.value.to_u64_digits();
    digits.resize(4, 0);
    digits.iter().flat_map(|d| d.to_le_bytes()).collect()
}

fn fields_to_bytes(xs: &Vec<F>) -> Vec<u8> {
    xs.iter().flat_map(|x| field_to_bytes(x)).collect()
}

fn fields_from_bytes(bytes: &[u8]) -> Vec<F> {
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

fn evals_over_boolean_hypercube(f: &MPolynomial<F>, result: &mut Vec<F>) {
    let el = f.get_num_vars();
    let comb = BitCombinationsDictOrder::new(el);
    for c in comb {
        let c_casted = c.iter().map(|v| F::from_value(*v)).collect::<Vec<_>>();
        result.push(f.evaluate(&c_casted));
    }
}

fn fold_factors_pointwise_cpu(evals: &mut [F], domain_size: usize, num_factors: usize) {
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
fn fold_into_half_cpu(
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
fn eval_folded_poly_cpu(
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

fn sum(data: &[F], len: usize) -> F {
    let mut result = F::zero();
    for i in 0..len {
        result = result.add_ref(&data[i]);
    }
    result
}

/// A CPU-based implementation of the Sum-Check protocol prover.
pub struct SumCheckProverCPU;

impl SumCheckProverCPU {
    /// Creates a new CPU-based Sum-Check prover.
    pub fn new() -> Self {
        Self
    }

    /// Computes the evaluations of the round polynomial `s_i(c)` for `c` in `{0, 1, ..., d}`.
    fn compute_round_polynomial_evals_cpu(
        &self,
        max_degree: usize,
        num_factors: usize,
        evals: &[F],
        num_remaining_vars: usize,
        domain_size: usize,
    ) -> Result<Vec<F>, Box<dyn std::error::Error>> {
        let s_eval_points: Vec<F> = (0..=max_degree).map(|d| F::from_value(d)).collect();
        let mut s_i_evaluations = Vec::with_capacity(max_degree + 1);

        for eval_point in s_eval_points.iter() {
            // Step 1: Evaluate the polynomial at `eval_point` for the current variable.
            let mut folded_evals = eval_folded_poly_cpu(
                evals,
                num_factors,
                domain_size,
                num_remaining_vars,
                eval_point,
            );

            // Step 2: Fold the factors of the newly evaluated polynomial.
            let current_domain_size = 1 << (num_remaining_vars - 1);
            fold_factors_pointwise_cpu(&mut folded_evals, current_domain_size, num_factors);

            // Step 3: Sum the results to get the evaluation s_i(eval_point).
            let round_sum: F = sum(&folded_evals, current_domain_size); //folded_evals.iter().take(current_domain_size).sum();
            s_i_evaluations.push(round_sum);
        }

        Ok(s_i_evaluations)
    }

    /// Executes the Sum-Check proving algorithm.
    pub fn prove(
        &self,
        max_degree: usize,
        polynomial_factors: &[MPolynomial<F>],
    ) -> Result<(F, Vec<u8>), Box<dyn std::error::Error>> {
        let num_variables = polynomial_factors
            .iter()
            .map(|p| p.get_num_vars())
            .max()
            .unwrap_or(0);
        let num_factors = polynomial_factors.len();

        // Generate the full evaluation table for all factors over the boolean hypercube.
        let mut evaluation_table = vec![];
        for f in polynomial_factors {
            evals_over_boolean_hypercube(f, &mut evaluation_table);
        }

        // Initialize transcript and add public parameters.
        let mut transcript = FiatShamirTransformer::new();
        transcript.push(&vec![bincode::serialize(&max_degree)?]);
        transcript.push(&vec![bincode::serialize(&num_factors)?]);
        transcript.push(&vec![bincode::serialize(&num_variables)?]);
        for factor in polynomial_factors {
            transcript.push(&vec![bincode::serialize(factor)?]);
        }

        let mut round_polynomials = Vec::with_capacity(num_variables);
        let mut challenges = Vec::with_capacity(num_variables);
        let mut num_remaining_vars = num_variables;
        let domain_size = 1 << num_variables;
        let s_eval_points: Vec<F> = (0..=max_degree).map(|d| F::from_value(d)).collect();

        // --- Calculate Initial Sum ---
        let claimed_sum = {
            let mut temp_evals = evaluation_table.clone();
            fold_factors_pointwise_cpu(&mut temp_evals, domain_size, num_factors);
            sum(&temp_evals, domain_size)
        };

        // --- Start Proving Rounds ---
        let mut current_evals = evaluation_table;

        for _round in 0..num_variables {
            // Step 1: Compute the round polynomial s_i's evaluations.
            let s_i_evaluations = self.compute_round_polynomial_evals_cpu(
                max_degree,
                num_factors,
                &current_evals,
                num_remaining_vars,
                domain_size,
            )?;

            // Step 2: Add evaluations to the transcript and interpolate the polynomial.
            for eval in &s_i_evaluations {
                transcript.push(&vec![bincode::serialize(eval)?]);
            }
            let s_i_poly = Polynomial::interpolate(&s_eval_points, &s_i_evaluations);
            round_polynomials.push(s_i_poly);

            // Step 3: Get the verifier's challenge for this round.
            let challenge = F::sample(&transcript.prover_fiat_shamir(32));
            challenges.push(challenge.clone());

            // Step 4: Fold the evaluation tables using the challenge for the next round.
            fold_into_half_cpu(
                &mut current_evals,
                num_factors,
                domain_size,
                num_remaining_vars,
                challenge,
            );

            num_remaining_vars -= 1;
        }

        Ok((claimed_sum, transcript.serialize()))
    }
}

struct SumCheckVerifier;

impl SumCheckVerifier {
    pub fn verify(
        max_degree: usize,
        polynomial_factors: &[MPolynomial<F>],
        claimed_sum: F,
        proof: &Vec<u8>,
    ) -> Result<bool, Box<dyn std::error::Error>> {
        let mut transcript = FiatShamirTransformer::deserialize(&proof);

        let num_variables = polynomial_factors
            .iter()
            .map(|p| p.get_num_vars())
            .max()
            .unwrap_or(0);
        let polynomial_product = polynomial_factors
            .iter()
            .skip(1)
            .fold(polynomial_factors[0].clone(), |acc, p| &acc * p);
        let s_eval_point: Vec<F> = (0..(max_degree + 1)).map(|d| F::from_value(d)).collect();

        let recovered_max_degree = transcript.pull();
        assert_eq!(
            vec![bincode::serialize(&max_degree).expect("Serialization failed")],
            recovered_max_degree
        );
        let recovered_num_factors = bincode::deserialize::<usize>(&transcript.pull()[0])?;
        assert_eq!(recovered_num_factors, polynomial_factors.len());
        let recovered_num_variables = bincode::deserialize::<usize>(&transcript.pull()[0])?;
        assert_eq!(recovered_num_variables, num_variables);
        for i in 0..recovered_num_factors {
            assert_eq!(
                vec![bincode::serialize(&polynomial_factors[i]).expect("Serialization failed")],
                transcript.pull()
            );
        }

        let mut s_polys = vec![];
        let mut challenges: Vec<F> = vec![];
        for i in 0..num_variables {
            let mut s_evals = vec![];
            for d in 0..(max_degree + 1) {
                let tmp_bytes = transcript.pull();
                let se = bincode::deserialize::<F>(&tmp_bytes[0])?;
                s_evals.push(se);
            }

            let s_poly = Polynomial::interpolate(&s_eval_point, &s_evals);
            s_polys.push(s_poly);

            let challenge = F::sample(&transcript.verifier_fiat_shamir(32));
            challenges.push(challenge);

            if i == 0 {
                assert_eq!(s_evals[0].add_ref(&s_evals[1]), claimed_sum);
            } else {
                assert_eq!(
                    s_evals[0].add_ref(&s_evals[1]),
                    s_polys[i - 1].eval(&challenges[i - 1]).sanitize()
                );
            }
        }
        assert_eq!(
            polynomial_product.evaluate(&challenges),
            s_polys[num_variables - 1]
                .eval(&challenges[num_variables - 1])
                .sanitize()
        );

        Ok(true)
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Starting Sumcheck Protocol Example (CPU)...");
    println!("-------------------------------------------");

    let mut dict_1 = HashMap::new();
    dict_1.insert(vec![0, 0, 0], F::from_value(1));
    dict_1.insert(vec![1, 0, 0], F::from_value(2));
    let factor_1 = MPolynomial::new(dict_1);

    let mut dict_2 = HashMap::new();
    dict_2.insert(vec![0, 0, 0], F::from_value(2));
    dict_2.insert(vec![0, 1, 0], F::from_value(3));
    let factor_2 = MPolynomial::new(dict_2);

    let mut dict_3 = HashMap::new();
    dict_3.insert(vec![0, 0, 0], F::from_value(3));
    dict_3.insert(vec![0, 0, 1], F::from_value(4));
    let factor_3 = MPolynomial::new(dict_3);

    let factors = vec![factor_1, factor_2, factor_3];

    let max_degree = 3;

    let prover = SumCheckProverCPU::new();
    println!("    Prover created. Generating proof...");
    let (claimed_sum, mut proof) = prover.prove(max_degree, &factors)?;
    println!("    ✅ Proof generated!");
    println!("    Prover's Claimed Sum (C1): {}", claimed_sum);

    debug_assert_eq!(
        sum_over_boolean_hypercube(
            &factors
                .iter()
                .skip(1)
                .fold(factors[0].clone(), |acc, p| &acc * p)
        ),
        claimed_sum
    );

    println!("    Verifier is now checking the proof interactively...");
    let is_valid = SumCheckVerifier::verify(max_degree, &factors, claimed_sum, &proof)?;
    println!("    ✅ Verification complete.");

    if is_valid {
        println!("\n🎉 SUCCESS: Proof is valid! 🎉");
    } else {
        // This branch will cause a panic due to the assert! above
        println!("\n❌ FAILURE: Proof is invalid! ❌");
    }

    Ok(())
}
