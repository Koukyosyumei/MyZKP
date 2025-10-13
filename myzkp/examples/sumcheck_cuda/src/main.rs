use std::fs;
use std::hash::Hash;
use std::collections::HashMap;

use cudarc;
use cudarc::nvrtc::Ptx;
use cudarc::nvrtc::compile_ptx;
use cudarc::driver::PushKernelArg;
use cudarc::driver::LaunchConfig;
use num_traits::identities::One;
use num_bigint::{BigInt, Sign};

use myzkp::modules::algebra::field::{FiniteFieldElement, ModEIP197};
use myzkp::modules::algebra::ring::Ring;
use myzkp::modules::algebra::mpolynomials::MPolynomial;
use myzkp::modules::algebra::polynomial::Polynomial;
use myzkp::modules::algebra::sumcheck::{BitCombinations, sum_over_boolean_hypercube};

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
    const FIELD_BYTES: usize = 256;

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
    let comb = BitCombinations::new(el, 0);
    for c in comb {
        let c_casted = c.iter().map(|v| F::from_value(*v)).collect::<Vec<_>>();
        result.push(f.evaluate(&c_casted));
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // let ptx = compile_ptx("").unwrap();
    let ptx = Ptx::from_file("../../src/modules/algebra/cuda/kernels/sumcheck.ptx"); //.unwrap();

    let a = F::one();
    let b = F::one();
    let x = vec![a, b];

    let ctx = cudarc::driver::CudaContext::new(0)?;
    let stream = ctx.default_stream();

    let module = ctx.load_module(ptx)?;
    let sum_kernel = module.load_function("sum")?;

    
    // copy a rust slice to the device
    let x_bytes = vec_f_to_bytes(&x);
    println!("x :{:?}", x_bytes);
    let inp = stream.memcpy_stod(&x_bytes)?;
    let mut out = stream.alloc_zeros::<u8>(32)?;
    

    let mut builder = stream.launch_builder(&sum_kernel);
    builder.arg(&inp);
    builder.arg(&mut out);
    builder.arg(&2usize);
    builder.arg(&0usize);
    unsafe { builder.launch(LaunchConfig::for_num_elems(2)) }?;
    let out_host: Vec<F> = vec_f_from_bytes(&stream.memcpy_dtov(&out)?);
    println!("out_host: {:?}", out_host);

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
    println!("g: {}", g);
    let h = sum_over_boolean_hypercube(&g);

    let mut evals = vec![];
    evals_over_boolean_hypercube(&factor_1, &mut evals);
    evals_over_boolean_hypercube(&factor_2, &mut evals);
    evals_over_boolean_hypercube(&factor_3, &mut evals);

    /*
    // or allocate directly

    let ptx = cudarc::nvrtc::compile_ptx("
    extern \"C\" __global__ void sin_kernel(float *out, const float *inp, const size_t numel) {
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < numel) {
        out[i] = sin(inp[i]);
    }
    }")?;

    // Dynamically load it into the device
    let module = ctx.load_module(ptx)?;
    let sin_kernel = module.load_function("sin_kernel")?;

    let mut builder = stream.launch_builder(&sin_kernel);
    builder.arg(&mut out);
    builder.arg(&inp);
    builder.arg(&100usize);
    unsafe { builder.launch(LaunchConfig::for_num_elems(100)) }?;

    let out_host: Vec<f32> = stream.memcpy_dtov(&out)?;
    assert_eq!(out_host, [1.0; 100].map(f32::sin));
    */

    Ok(())
}
