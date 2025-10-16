use std::collections::HashMap;
use std::sync::Arc;
use std::time;

use cudarc;
use cudarc::driver::CudaStream;
use cudarc::driver::LaunchConfig;
use cudarc::driver::{PushKernelArg, CudaFunction, CudaSlice};
use cudarc::nvrtc::Ptx;
use num_bigint::{BigInt, Sign};
use rand::Rng;
use rand::{rngs::StdRng, SeedableRng};

use myzkp::modules::algebra::fiat_shamir::FiatShamirTransformer;
use myzkp::modules::algebra::field::{Field, FiniteFieldElement, ModEIP197};
use myzkp::modules::algebra::mpolynomials::MPolynomial;
use myzkp::modules::algebra::polynomial::Polynomial;
use myzkp::modules::algebra::ring::Ring;
use myzkp::modules::algebra::sumcheck::sum_over_boolean_hypercube;

use sumcheck_cuda::prover::{CudaBackend, SumCheckProverCPU, SumCheckProverGPU};
use sumcheck_cuda::verifier::SumCheckVerifier;
use sumcheck_cuda::utils::{F, BitCombinationsDictOrder};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Starting Sumcheck Protocol Example...");
    println!("-------------------------------------------");
    let is_gpu = true;

    let num_vars = 11;
    let num_factors = 2;
    let max_degree = num_factors;

    let mut factors = vec![];
    let mut rng = StdRng::seed_from_u64(42);
    for _ in 0..num_factors {
        let mut dict = HashMap::new();
        let comb = BitCombinationsDictOrder::new(num_vars);
        for c in comb {
            let c_casted: Vec<usize> = c.iter().map(|&b| b as usize).collect();
            dict.insert(c_casted, F::from_value(rng.random_range(0..256)));
        }
        factors.push(MPolynomial::new(dict));
    }

    let num_blocks_per_poly = 1;
    let num_threads_per_block = 256;

    println!("    Prover created. Generating proof...");
    let start_time = time::Instant::now();    
    let (claimed_sum, mut proof) = {
        if is_gpu {    
            println!("    Initializing CUDA Backend...");
            let cuda_backend = CudaBackend::new(0)?;
            println!("    ✅ CUDA Backend initialized.");
            let prover = SumCheckProverGPU::new(
                num_blocks_per_poly,
                num_threads_per_block,
                &cuda_backend,
            );
            prover.prove(max_degree, &factors)?
        } else {
            let prover = SumCheckProverCPU::new();
            prover.prove(max_degree, &factors)?
        }
    };
    println!("    ✅ Proof generated!");
    println!("        - Prover's Claimed Sum (C1): {}", claimed_sum);
    println!("        - Proving Time: {:?}", start_time.elapsed());

    debug_assert_eq!(sum_over_boolean_hypercube(&factors.iter().skip(1).fold(factors[0].clone(), |acc, p| &acc * p)), claimed_sum);
    
    println!("    Verifier is now checking the proof interactively...");
    let start_time = time::Instant::now();
    let is_valid = SumCheckVerifier::verify(max_degree, &factors, claimed_sum, &proof)?;
    println!("    ✅ Verification complete.");
    println!("        - Verification Time: {:?}", start_time.elapsed());

    if is_valid {
        println!("\n🎉 SUCCESS: Proof is valid! 🎉");
    } else {
        // This branch will cause a panic due to the assert! above
        println!("\n❌ FAILURE: Proof is invalid! ❌");
    }

    Ok(())
}
