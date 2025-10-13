use std::fs;
use std::hash::Hash;
use std::collections::HashMap;

use cudarc;
use cudarc::nvrtc::Ptx;
use cudarc::nvrtc::compile_ptx;
use cudarc::driver::PushKernelArg;
use cudarc::driver::LaunchConfig;
use num_traits::identities::One;
use num_traits::Zero;
use num_bigint::{BigInt, Sign};

use myzkp::modules::algebra::field::{FiniteFieldElement, ModEIP197};
use myzkp::modules::algebra::ring::Ring;
use myzkp::modules::algebra::field::Field;
use myzkp::modules::algebra::mpolynomials::MPolynomial;
use myzkp::modules::algebra::polynomial::Polynomial;
use myzkp::modules::algebra::sumcheck::{sum_over_boolean_hypercube};
use myzkp::modules::algebra::fiat_shamir::FiatShamirTransformer;

#[derive(Debug)]
struct BitCombinations {
    len: usize,
    current: usize,
}

impl BitCombinations {
    pub fn new(length: usize) -> Self {
        assert!(
            length <= usize::BITS as usize,
            "Length exceeds available bits"
        );
        BitCombinations { len: length, current: 0 }
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
        for i in (0..self.len).rev() {
            let bit = (n >> i) & 1;
            combination.push(bit as u8);
        }

        self.current += 1;
        Some(combination)
    }
}

type F = FiniteFieldElement<ModEIP197>;

fn f_to_bytes(x: &F) -> Vec<u8> {
    let (_, mut digits) = x.value.to_u64_digits();
    digits.resize(4, 0);
    digits.iter().flat_map(|d| d.to_le_bytes()).collect()
}

fn vec_f_to_bytes(xs: &Vec<F>) -> Vec<u8> {
    xs.iter().flat_map(|x| f_to_bytes(x)).collect()
}

fn vec_f_from_bytes(bytes: &[u8]) -> Vec<F> {
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
    let comb = BitCombinations::new(el);
    for c in comb {
        let c_casted = c.iter().map(|v| F::from_value(*v)).collect::<Vec<_>>();
        result.push(f.evaluate(&c_casted));
    }
}

//fn prove_sumcheck() {
//
//}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let ptx = Ptx::from_file("../../src/modules/algebra/cuda/kernels/sumcheck.ptx"); //.unwrap();

    let ctx = cudarc::driver::CudaContext::new(0)?;
    let stream = ctx.default_stream();

    let module = ctx.load_module(ptx)?;
    let fold_factors_pointwise_kernel = module.load_function("fold_factors_pointwise")?;
    let fold_into_half_kernel = module.load_function("fold_into_half")?;
    let eval_folded_poly_kernel = module.load_function("eval_folded_poly")?;
    let sum_kernel = module.load_function("sum")?;

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

    let g = &(&factor_1 * &factor_2) * &factor_3;
    println!("f_1: {}", factor_1);
    println!("f_2: {}", factor_2);
    println!("f_3: {}", factor_3);
    println!("g: {}", g);
    let h = sum_over_boolean_hypercube(&g);
    println!("h: {}", h);

    let mut evals = vec![];
    evals_over_boolean_hypercube(&factor_1, &mut evals);
    evals_over_boolean_hypercube(&factor_2, &mut evals);
    evals_over_boolean_hypercube(&factor_3, &mut evals);
    let evals_bytes = vec_f_to_bytes(&evals);
    for i in 0..3 {
        for j in 0..8 {
            print!("{}, ", evals[i * 8 + j]);
        }
        print!("\n");
    }

    let max_degree = 3;
    let num_factors: usize = 3;
    let num_vars = 3;
    let domain_size = 1 << num_vars;
    
    let s_eval_point: Vec<F> = (0..(max_degree+1)).map(|d| F::from_value(d)).collect();
    let mut num_remaining_vars = num_vars;
    let mut s_polys = vec![];
    let mut challenges: Vec<F> = vec![];
    
    let num_blocks_per_poly = 1;
    let num_threads_per_block = 256;

    let launch_config = LaunchConfig {
        grid_dim: (num_blocks_per_poly * (num_factors as u32), 1, 1),
        block_dim: (num_threads_per_block, 1, 1),
        shared_mem_bytes: 128,
    };

    let mut proof_stream = FiatShamirTransformer::new();
    proof_stream.push(&vec![bincode::serialize(&max_degree).expect("Serialization failed")]);
    proof_stream.push(&vec![bincode::serialize(&num_vars).expect("Serialization failed")]);
    proof_stream.push(&vec![bincode::serialize(&g).expect("Serialization failed")]);

    let mut evals_dev = stream.memcpy_stod(&evals_bytes)?;
    let mut s_evals_dev = stream.alloc_zeros::<u8>(32 * (max_degree + 1))?;
    let mut buf_dev = stream.alloc_zeros::<u8>(32 * ((1 << (num_remaining_vars - 1)) * num_factors))?;

    for i in 0..num_vars {
        for d in 0..(max_degree+1) {
            let eval_point = F::from_value(d);
            let eval_point_bytes = f_to_bytes(&eval_point);
            let eval_point_dev = stream.memcpy_stod(&eval_point_bytes)?;
        
            let mut builder = stream.launch_builder(&eval_folded_poly_kernel);
            builder.arg(&num_remaining_vars);
            builder.arg(&domain_size);
            builder.arg(&num_blocks_per_poly);
            builder.arg(&evals_dev);
            builder.arg(&mut buf_dev);
            builder.arg(&eval_point_dev);
            unsafe { builder.launch(launch_config) }?;

            let mut builder = stream.launch_builder(&fold_factors_pointwise_kernel);
            builder.arg(&mut buf_dev);
            let half_domain_size: usize = 1 << (num_remaining_vars - 1);
            builder.arg(&half_domain_size);
            builder.arg(&num_factors);
            unsafe { builder.launch(launch_config) }?;

            let half_half_domain_size: usize = half_domain_size >> 1;
            let mut builder = stream.launch_builder(&sum_kernel);
            builder.arg(&buf_dev);
            builder.arg(&mut s_evals_dev);
            builder.arg(&half_half_domain_size);
            builder.arg(&d);
            unsafe { builder.launch(launch_config) }?;
        }

        let s_evals_host: Vec<F> = vec_f_from_bytes(&stream.memcpy_dtov(&s_evals_dev)?);
        for se in &s_evals_host {
            proof_stream.push(&vec![bincode::serialize(&se).expect("Serialization failed")]);
        }

        let s_poly = Polynomial::interpolate(&s_eval_point, &s_evals_host);
        s_polys.push(s_poly);

        let challenge = F::sample(&proof_stream.prover_fiat_shamir(32));
        let challenge_bytes = f_to_bytes(&challenge);
        let challenge_dev = stream.memcpy_stod(&challenge_bytes)?;
        let mut builder = stream.launch_builder(&fold_into_half_kernel);
        builder.arg(&num_remaining_vars);
        builder.arg(&domain_size);
        builder.arg(&num_blocks_per_poly);
        builder.arg(&mut evals_dev);
        builder.arg(&challenge_dev);
        unsafe { builder.launch(launch_config) }?;
        challenges.push(challenge);

        if i == 0 {
            assert_eq!(s_evals_host[0].add_ref(&s_evals_host[1]), h);
        } else {
            assert_eq!(s_evals_host[0].add_ref(&s_evals_host[1]), s_polys[i-1].eval(&challenges[i-1]).sanitize());
        }
        num_remaining_vars -= 1;
    }

    assert_eq!(g.evaluate(&challenges), s_polys[num_vars - 1].eval(&challenges[num_vars - 1]).sanitize());

    let recovered_max_degree = proof_stream.pull();   
    assert_eq!(vec![bincode::serialize(&max_degree).expect("Serialization failed")], recovered_max_degree);
    let recovered_num_vars = proof_stream.pull();
    assert_eq!(vec![bincode::serialize(&num_vars).expect("Serialization failed")], recovered_num_vars);
    let recovered_g = proof_stream.pull();
    assert_eq!(vec![bincode::serialize(&g).expect("Serialization failed")], recovered_g);

    let mut s_polys = vec![];
    let mut challenges: Vec<F> = vec![];
    for i in 0..num_vars {
        let mut s_evals = vec![];
        for d in 0..(max_degree+1) {
            let tmp_bytes = proof_stream.pull();
            let se = bincode::deserialize::<F>(&tmp_bytes[0])?;
            s_evals.push(se);
        }

        let s_poly = Polynomial::interpolate(&s_eval_point, &s_evals);
        s_polys.push(s_poly);

        let challenge = F::sample(&proof_stream.verifier_fiat_shamir(32));
        challenges.push(challenge);

        if i == 0 {
            assert_eq!(s_evals[0].add_ref(&s_evals[1]), h);
        } else {
            assert_eq!(s_evals[0].add_ref(&s_evals[1]), s_polys[i-1].eval(&challenges[i-1]).sanitize());
        }
    }


    Ok(())
}
