use std::sync::Arc;

use cudarc;
use cudarc::driver::CudaStream;
use cudarc::driver::LaunchConfig;
use cudarc::driver::{PushKernelArg, CudaFunction, CudaSlice};
use cudarc::nvrtc::Ptx;

use myzkp::modules::algebra::fiat_shamir::FiatShamirTransformer;
use myzkp::modules::algebra::field::Field;
use myzkp::modules::algebra::mpolynomials::MPolynomial;
use myzkp::modules::algebra::polynomial::Polynomial;
use myzkp::modules::algebra::ring::Ring;

use crate::utils::{fold_factors_pointwise_cpu, fold_into_half_cpu, eval_folded_poly_cpu, sum_cpu, field_to_bytes, fields_to_bytes, fields_from_bytes, F};

pub struct CudaBackend {
    stream: Arc<CudaStream>,
    eval_all_binary_combinations_kernel: CudaFunction,
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
            eval_all_binary_combinations_kernel: module.load_function("eval_all_binary_combinations")?,
            fold_factors_pointwise_kernel: module.load_function("fold_factors_pointwise")?,
            fold_into_half_kernel: module.load_function("fold_into_half")?,
            eval_folded_poly_kernel: module.load_function("eval_folded_poly")?,
            sum_kernel: module.load_function("sum")?,
        })
    }
}

pub struct MPolyKernelInput {
    pub coeffs: Vec<F>,
    pub expo_flat: Vec<u32>,
    pub offsets: Vec<u32>,
    pub lens: Vec<u32>,
}

pub fn convert_mpoly_to_kernel_format(mpoly: &MPolynomial<F>) -> MPolyKernelInput {
    let mut coeffs = Vec::new();
    let mut expo_flat = Vec::new();
    let mut offsets = Vec::new();
    let mut lens = Vec::new();
    let mut current_offset = 0u32;

    for (exponents, coeff) in &mpoly.dictionary {
        coeffs.push(coeff.clone());
        offsets.push(current_offset);
        lens.push(exponents.len() as u32);

        for &e in exponents {
            expo_flat.push(e as u32);
        }

        current_offset += exponents.len() as u32;
    }

    MPolyKernelInput {
        coeffs,
        expo_flat,
        offsets,
        lens,
    }
}

pub struct SumCheckProverGPU<'a> {
    num_blocks_per_poly: usize,
    num_threads_per_block: usize,
    cuda_backend: &'a CudaBackend,
}

impl<'a> SumCheckProverGPU<'a> {
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
        /*
        let mut evaluation_table = vec![];
        for f in polynomial_factors {
            evals_over_boolean_hypercube(&f, &mut evaluation_table);
        }
        */

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


        // Evaluation Table
        let mut evaluation_table_dev_p = self
            .cuda_backend
            .stream
            .alloc_zeros::<u8>(32 * domain_size * num_factors)?;
        for i in 0..num_factors {
            let mpoly_info = convert_mpoly_to_kernel_format(&polynomial_factors[i]);
            let coeffs_bytes = fields_to_bytes(&mpoly_info.coeffs);
            let coeffs_bytes_dev = self.cuda_backend.stream.memcpy_stod(&coeffs_bytes)?;
            let expo_flat_dev = self.cuda_backend.stream.memcpy_stod(&mpoly_info.expo_flat)?;
            let offsets_dev = self.cuda_backend.stream.memcpy_stod(&mpoly_info.offsets)?;
            let lens_dev = self.cuda_backend.stream.memcpy_stod(&mpoly_info.lens)?;
            let offset_p = (i * domain_size) as usize;        

            let mut builder = self
                .cuda_backend
                .stream
                .launch_builder(&self.cuda_backend.eval_all_binary_combinations_kernel);
            builder.arg(&mut evaluation_table_dev_p);
            builder.arg(&offset_p);
            builder.arg(&num_variables);
            builder.arg(&domain_size);
            builder.arg(&coeffs_bytes_dev);
            builder.arg(&expo_flat_dev);
            builder.arg(&offsets_dev);
            builder.arg(&lens_dev);
            unsafe { builder.launch(launch_config) }?;
        }
        let evaluation_table: Vec<F> =
            fields_from_bytes(&self.cuda_backend.stream.memcpy_dtov(&evaluation_table_dev_p)?);
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
        builder.arg(&mut sum_result_dev);
        builder.arg(&tmp_evals_dev);
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
            builder.arg(&mut buf_dev);
            builder.arg(evals_dev);
            builder.arg(&eval_point_dev);
            builder.arg(&num_remaining_vars);
            builder.arg(&domain_size);
            builder.arg(&self.num_blocks_per_poly);
            unsafe { builder.launch(*launch_config) }?;

            let mut builder = self
                .cuda_backend
                .stream
                .launch_builder(&self.cuda_backend.fold_factors_pointwise_kernel);
            let half_domain_size: usize = 1 << (num_remaining_vars - 1);
            builder.arg(&mut buf_dev);
            builder.arg(&half_domain_size);
            builder.arg(&num_factors);
            unsafe { builder.launch(*launch_config) }?;

            let half_half_domain_size: usize = half_domain_size >> 1;
            let mut builder = self
                .cuda_backend
                .stream
                .launch_builder(&self.cuda_backend.sum_kernel);
            builder.arg(&mut *s_evals_dev);
            builder.arg(&buf_dev);
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
        builder.arg(evals_dev);
        builder.arg(&domain_size);
        builder.arg(&num_remaining_vars);
        builder.arg(&challenge_dev);
        builder.arg(&self.num_blocks_per_poly);
        unsafe { builder.launch(*launch_config) }?;
        Ok(())
    }
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
            let round_sum: F = sum_cpu(&folded_evals, current_domain_size); //folded_evals.iter().take(current_domain_size).sum();
            s_i_evaluations.push(round_sum);
        }

        Ok(s_i_evaluations)
    }

    /// Executes the Sum-Check proving algorithm.
    pub fn prove(
        &self,
        evaluation_table: &mut Vec<F>,
        max_degree: usize,
        polynomial_factors: &[MPolynomial<F>]
    ) -> Result<(F, Vec<u8>), Box<dyn std::error::Error>> {
        let num_variables = polynomial_factors
            .iter()
            .map(|p| p.get_num_vars())
            .max()
            .unwrap_or(0);
        let num_factors = polynomial_factors.len();

        /*
        // Generate the full evaluation table for all factors over the boolean hypercube.
        let mut evaluation_table = vec![];
        for f in polynomial_factors {
            evals_over_boolean_hypercube(f, &mut evaluation_table);
        }
        */

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
            sum_cpu(&temp_evals, domain_size)
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
