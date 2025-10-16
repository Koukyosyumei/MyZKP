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

use sumcheck_cuda::prover::{CudaBackend, SumCheckProverCPU, SumCheckProverGPU};
use sumcheck_cuda::verifier::SumCheckVerifier;
use sumcheck_cuda::utils::{F};

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

    let prover = SumCheckProverGPU::new(
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
