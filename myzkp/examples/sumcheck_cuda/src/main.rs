use std::collections::HashMap;
use std::sync::Arc;

use cudarc;
use cudarc::driver::CudaStream;
use cudarc::driver::LaunchConfig;
use cudarc::driver::{PushKernelArg, CudaFunction, CudaSlice};
use cudarc::nvrtc::Ptx;
use num_bigint::{BigInt, Sign};

use myzkp::modules::algebra::fiat_shamir::FiatShamirTransformer;
use myzkp::modules::algebra::field::{Field, FiniteFieldElement, ModEIP197};
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

struct CudaBackend {
    stream: Arc<CudaStream>,
    fold_factors_pointwise_kernel: CudaFunction,
    fold_into_half_kernel: CudaFunction,
    eval_folded_poly_kernel: CudaFunction,
    sum_kernel: CudaFunction,
}

impl CudaBackend {
    pub fn new(device_id: usize) -> Result<Self, Box<dyn std::error::Error>> {
        let ctx = cudarc::driver::CudaContext::new(device_id)?;
        let stream = ctx.default_stream();
        let ptx = Ptx::from_file("../../src/modules/algebra/cuda/kernels/sumcheck.ptx");
        let module = ctx.load_module(ptx)?;

        Ok(Self {
            stream,
            fold_factors_pointwise_kernel: module.load_function("fold_factors_pointwise")?,
            fold_into_half_kernel: module.load_function("fold_into_half")?,
            eval_folded_poly_kernel: module.load_function("eval_folded_poly")?,
            sum_kernel: module.load_function("sum")?,
        })
    }
}

struct SumCheckProver<'a> {
    num_blocks_per_poly: usize,
    num_threads_per_block: usize,
    cuda_backend: &'a CudaBackend,
}

impl<'a> SumCheckProver<'a> {
    pub fn new(
        num_blocks_per_poly: usize,
        num_threads_per_block: usize,
        cuda_backend: &'a CudaBackend,
    ) -> Self {
        Self {
            num_blocks_per_poly,
            num_threads_per_block,
            cuda_backend,
        }
    }

    pub fn prove(&self, max_degree: usize, polynomial_factors: &[MPolynomial<F>]) -> Result<(F, Vec<u8>), Box<dyn std::error::Error>> {
        let num_variables = polynomial_factors.iter().map(|p| p.get_num_vars()).max().unwrap_or(0);
        let num_factors = polynomial_factors.len();
        let mut evaluation_table = vec![];
        for f in polynomial_factors {
            evals_over_boolean_hypercube(&f, &mut evaluation_table);
        }

        let mut transcript = FiatShamirTransformer::new();
        transcript.push(&vec![
            bincode::serialize(&max_degree).expect("Serialization failed")
        ]);
        transcript.push(&vec![
            bincode::serialize(&num_factors).expect("Serialization failed")
        ]);
        transcript.push(&vec![
            bincode::serialize(&num_variables).expect("Serialization failed")
        ]);
        for i in 0..num_factors {
            transcript.push(&vec![
                bincode::serialize(&polynomial_factors[i]).expect("Serialization failed")
            ]);
        }

        let mut round_polynomials = Vec::with_capacity(num_variables);
        let mut challenges = Vec::with_capacity(num_variables);
        let mut num_remaining_vars = num_variables;
        let domain_size = 1 << num_variables;
        let s_eval_points: Vec<F> = (0..=max_degree).map(|d| F::from_value(d)).collect();

        let launch_config = LaunchConfig {
            grid_dim: (
                (self.num_blocks_per_poly as u32) * (num_factors as u32),
                1,
                1,
            ),
            block_dim: (self.num_threads_per_block as u32, 1, 1),
            shared_mem_bytes: 128,
        };

        let all_evals_bytes = fields_to_bytes(&evaluation_table);

        // Calculate Sum
        let mut tmp_evals_dev = self.cuda_backend.stream.memcpy_stod(&all_evals_bytes)?;
        let mut sum_result_dev = self
            .cuda_backend
            .stream
            .alloc_zeros::<u8>(32)?;

        let mut builder = self
            .cuda_backend
            .stream
            .launch_builder(&self.cuda_backend.fold_factors_pointwise_kernel);
        builder.arg(&mut tmp_evals_dev);
        builder.arg(&domain_size);
        builder.arg(&num_factors);
        unsafe { builder.launch(launch_config) }?;

        let half_domain_size: usize = domain_size >> 1;
        let mut builder = self
            .cuda_backend
            .stream
            .launch_builder(&self.cuda_backend.sum_kernel);
        builder.arg(&tmp_evals_dev);
        builder.arg(&mut sum_result_dev);
        builder.arg(&half_domain_size);
        builder.arg(&0);
        unsafe { builder.launch(launch_config) }?;
        let sum_result_host: F =
            fields_from_bytes(&self.cuda_backend.stream.memcpy_dtov(&sum_result_dev)?)[0].clone();

        // Start Proving
        let mut evals_dev = self.cuda_backend.stream.memcpy_stod(&all_evals_bytes)?;
        let mut s_evals_dev = self
            .cuda_backend
            .stream
            .alloc_zeros::<u8>(32 * (max_degree + 1))?;

        for _round in 0..num_variables {
            // Step 1: Compute the round polynomial s_i.
            let s_i_evaluations = self.compute_round_polynomial_evals(
                max_degree,
                num_factors,
                &evals_dev,
                &mut s_evals_dev,
                num_remaining_vars,
                domain_size,
                &launch_config,
            )?;

            // Step 2: Add evaluations to the transcript.
            s_i_evaluations.iter().for_each(|e| {
                transcript.push(&vec![bincode::serialize(e).unwrap()]);
            });

            let s_i_poly = Polynomial::interpolate(&s_eval_points, &s_i_evaluations);
            round_polynomials.push(s_i_poly);

            // Step 3: Get the verifier's challenge for this round.
            let challenge = F::sample(&transcript.prover_fiat_shamir(32));
            challenges.push(challenge.clone());

            // Step 4: Fold the evaluation tables using the challenge.
            // This corresponds to the update rule A[x] <- (1-r)A[0,x] + rA[1,x][cite: 266].
            self.fold_evaluations(
                domain_size,
                &mut evals_dev,
                num_remaining_vars,
                challenge,
                &launch_config,
            )?;

            num_remaining_vars -= 1;
        }

        Ok((sum_result_host, transcript.serialize()))
    }

    /// Computes the evaluations of the round polynomial `s_i(c)` for `c` in `{0, 1, ..., d}`.
    fn compute_round_polynomial_evals(
        &self,
        max_degree: usize,
        num_factors: usize,
        evals_dev: &CudaSlice<u8>,
        s_evals_dev: &mut CudaSlice<u8>,
        num_remaining_vars: usize,
        domain_size: usize,
        launch_config: &LaunchConfig,
    ) -> Result<Vec<F>, Box<dyn std::error::Error>> {
        let s_eval_points: Vec<F> = (0..=max_degree).map(|d| F::from_value(d)).collect();

        let mut buf_dev = self
            .cuda_backend
            .stream
            .alloc_zeros::<u8>(32 * ((1 << (num_remaining_vars - 1)) * num_factors))?;

        for (idx, eval_point) in s_eval_points.iter().enumerate() {
            let eval_point_bytes = field_to_bytes(eval_point);
            let eval_point_dev = self.cuda_backend.stream.memcpy_stod(&eval_point_bytes)?;

            let mut builder = self
                .cuda_backend
                .stream
                .launch_builder(&self.cuda_backend.eval_folded_poly_kernel);
            builder.arg(&num_remaining_vars);
            builder.arg(&domain_size);
            builder.arg(&self.num_blocks_per_poly);
            builder.arg(evals_dev);
            builder.arg(&mut buf_dev);
            builder.arg(&eval_point_dev);
            unsafe { builder.launch(*launch_config) }?;

            let mut builder = self
                .cuda_backend
                .stream
                .launch_builder(&self.cuda_backend.fold_factors_pointwise_kernel);
            builder.arg(&mut buf_dev);
            let half_domain_size: usize = 1 << (num_remaining_vars - 1);
            builder.arg(&half_domain_size);
            builder.arg(&num_factors);
            unsafe { builder.launch(*launch_config) }?;

            let half_half_domain_size: usize = half_domain_size >> 1;
            let mut builder = self
                .cuda_backend
                .stream
                .launch_builder(&self.cuda_backend.sum_kernel);
            builder.arg(&buf_dev);
            builder.arg(&mut *s_evals_dev);
            builder.arg(&half_half_domain_size);
            builder.arg(&idx);
            unsafe { builder.launch(*launch_config) }?;
        }

        // Retrieve results from GPU.
        let s_evals_host: Vec<F> =
            fields_from_bytes(&self.cuda_backend.stream.memcpy_dtov(s_evals_dev)?);
        Ok(s_evals_host)
    }

    fn fold_evaluations(
        &self,
        domain_size: usize,
        evals_dev: &mut CudaSlice<u8>,
        num_remaining_vars: usize,
        challenge: F,
        launch_config: &LaunchConfig,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let challenge_bytes = field_to_bytes(&challenge);
        let challenge_dev = self.cuda_backend.stream.memcpy_stod(&challenge_bytes)?;

        let mut builder = self
            .cuda_backend
            .stream
            .launch_builder(&self.cuda_backend.fold_into_half_kernel);
        builder.arg(&num_remaining_vars);
        builder.arg(&domain_size);
        builder.arg(&self.num_blocks_per_poly);
        builder.arg(evals_dev);
        builder.arg(&challenge_dev);
        unsafe { builder.launch(*launch_config) }?;
        Ok(())
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

        let num_variables = polynomial_factors.iter().map(|p| p.get_num_vars()).max().unwrap_or(0);
        let polynomial_product = polynomial_factors.iter().skip(1).fold(polynomial_factors[0].clone(), |acc, p| &acc * p);
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
    println!("Starting Sumcheck Protocol Example...");
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
    let num_blocks_per_poly = 1;
    let num_threads_per_block = 256;

    println!("    Initializing CUDA Backend...");
    let cuda_backend = CudaBackend::new(0)?;
    println!("    ✅ CUDA Backend initialized.");

    let prover = SumCheckProver::new(
        num_blocks_per_poly,
        num_threads_per_block,
        &cuda_backend,
    );
    println!("    Prover created. Generating proof...");
    let (claimed_sum, mut proof) = prover.prove(max_degree, &factors)?;
    println!("    ✅ Proof generated!");
    println!("    Prover's Claimed Sum (C1): {}", claimed_sum);

    debug_assert_eq!(sum_over_boolean_hypercube(&factors.iter().skip(1).fold(factors[0].clone(), |acc, p| &acc * p)), claimed_sum);
    
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
