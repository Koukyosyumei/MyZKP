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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    //let contents = fs::read_to_string("../../src/modules/algebra/cuda/kernels/sumcheck.hpp")?;
    let ptx = Ptx::from_file("../../src/modules/algebra/cuda/kernels/sumcheck.cu"); //.unwrap();

    let a = F::one();

    let ctx = cudarc::driver::CudaContext::new(0)?;
    let stream = ctx.default_stream();

    // copy a rust slice to the device
    let inp = stream.memcpy_stod(&[1.0f32; 100])?;

    // or allocate directly
    let mut out = stream.alloc_zeros::<f32>(100)?;

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

    Ok(())
}
