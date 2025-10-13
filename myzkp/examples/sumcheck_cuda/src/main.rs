use std::fs;

use cudarc;
use cudarc::nvrtc::Ptx;
use cudarc::nvrtc::compile_ptx;
use cudarc::driver::PushKernelArg;
use cudarc::driver::LaunchConfig;
use num_traits::identities::One;

use myzkp::modules::algebra::field::{FiniteFieldElement, ModEIP197};

type F = FiniteFieldElement<ModEIP197>;

fn f_to_bytes(x: &F) -> Vec<u8> {
    let (_, digits) = x.value.to_u64_digits();
    digits.iter().flat_map(|d| d.to_le_bytes()).collect()
}

fn vec_f_to_bytes(xs: &Vec<F>) -> Vec<u8> {
    xs.iter().flat_map(|x| f_to_bytes(x)).collect()
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
    // println!(":{}")
    let inp = stream.memcpy_stod(&x_bytes)?;
    let mut out = stream.alloc_zeros::<u8>(32)?;
    

    let mut builder = stream.launch_builder(&sum_kernel);
    builder.arg(&inp);
    builder.arg(&mut out);
    builder.arg(&2usize);
    builder.arg(&0usize);
    unsafe { builder.launch(LaunchConfig::for_num_elems(2)) }?;
    let out_host: Vec<u8> = stream.memcpy_dtov(&out)?;
    println!("out_host: {:?}", out_host);

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
