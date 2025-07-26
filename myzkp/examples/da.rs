use std::{cmp::min, env};

use myzkp::modules::das::{
    avail::Avail,
    celestia::Celestia,
    eigenda::EigenDA,
    utils::{reset_metrics, DataAvailabilitySystem, SamplePosition, METRICS},
};

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} [eigenda | celestia | avail]", args[0]);
        std::process::exit(1);
    }
    let target = &args[1];

    let data_sizes: Vec<usize> = vec![16, 64, 256, 1024];
    let sqrt_data_sizes: Vec<usize> = vec![4, 8, 16, 32];

    for (data_size, sqrt_data_size) in data_sizes.iter().zip(sqrt_data_sizes.iter()) {
        let data: Vec<u8> = (0..*data_size).map(|i| (i % 256) as u8).collect();

        match target.as_str() {
            "eigenda" => {
                println!("# EigenDA");
                let num_operators = 8;
                let num_verification = 5;
                let expansion_factor = 4.0;
                let chunk_size =
                    ((*data_size as f64) * expansion_factor / num_operators as f64) as usize;

                let params = EigenDA::setup(chunk_size, expansion_factor, *data_size);
                let encoded = EigenDA::encode(&data, &params);
                let commit = EigenDA::commit(&encoded, &params);

                for i in 0..(num_verification) {
                    let position = SamplePosition {
                        row: 0,
                        col: i,
                        is_row: false,
                    };
                    let isvalid = EigenDA::verify(&position, &encoded, &commit, &params);
                    assert!(isvalid);
                }
                METRICS.with(|metrics| {
                    println!("{:#?}", *metrics.borrow());
                });
                reset_metrics();
            }
            "celestia" => {
                println!("# Celestia");
                let expansion_factor = 2;
                let base_num_sampling = 16;

                let params = Celestia::setup(*sqrt_data_size, expansion_factor as f64, *data_size);
                let encoded = Celestia::encode(&data, &params);
                let commit = Celestia::commit(&encoded, &params);
                for i in 0..min(
                    (sqrt_data_size * expansion_factor).pow(2),
                    base_num_sampling,
                ) {
                    let position = SamplePosition {
                        row: i / (sqrt_data_size * expansion_factor),
                        col: i % (sqrt_data_size * expansion_factor),
                        is_row: false,
                    };
                    assert!(Celestia::verify(&position, &encoded, &commit, &params));
                }
                METRICS.with(|metrics| {
                    println!("{:#?}", *metrics.borrow());
                });
                reset_metrics();
            }
            "avail" => {
                println!("# Avail");
                let expansion_factor = 2;
                let chunk_size = 8;
                let base_num_sampling = 8;

                let params = Avail::setup(chunk_size, expansion_factor as f64, *data_size);
                let encoded = Avail::encode(&data, &params);
                let commit = Avail::commit(&encoded, &params);
                for i in 0..min(chunk_size * expansion_factor as usize, base_num_sampling) {
                    let position0 = SamplePosition {
                        row: 0,
                        col: i,
                        is_row: false,
                    };
                    assert!(Avail::verify(&position0, &encoded, &commit, &params));
                }
                METRICS.with(|metrics| {
                    println!("{:#?}", *metrics.borrow());
                });
                reset_metrics();
            }
            _ => {
                eprintln!("Unknown target: {}", target);
                eprintln!("Usage: {} [eigenda | celestia | avail]", args[0]);
                std::process::exit(1);
            }
        }
    }
}
