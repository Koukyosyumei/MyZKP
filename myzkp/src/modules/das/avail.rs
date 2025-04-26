use std::time::Instant;

use crate::modules::algebra::curve::bn128::FqOrder;
use crate::modules::algebra::curve::bn128::BN128;
use crate::modules::algebra::kzg::{
    commit_kzg, open_kzg, setup_kzg, verify_kzg, CommitmentKZG, ProofKZG, PublicKeyKZG,
};
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::reedsolomon::{decode_rs2d, encode_rs2d, setup_rs2d};
use crate::modules::algebra::ring::Ring;
use crate::modules::das::utils::{DataAvailabilitySystem, SamplePosition, SystemMetrics, METRICS};

pub struct EncodedDataAvail {
    pub codewords: Vec<Vec<u8>>,
    pub data_size: usize, // original data size
}

pub struct CommitmentAvail {
    pub row_commitments: Vec<CommitmentKZG>,
    pub col_commitments: Vec<CommitmentKZG>,
}

pub struct PublicParamsAvail {
    pub codeword_size: usize, // the expanded side length
    pub pk: PublicKeyKZG,
    pub chunk_size: usize, // the side length when reshaping the original into a 2D array
}

pub struct Avail;

impl DataAvailabilitySystem for Avail {
    type EncodedData = EncodedDataAvail;
    type Commitment = CommitmentAvail;
    type PublicParams = PublicParamsAvail;

    fn setup(chunk_size: usize, expansion_factor: f64) -> Self::PublicParams {
        let g1 = BN128::generator_g1();
        let g2 = BN128::generator_g2();
        let codeword_size = (chunk_size as f64 * expansion_factor.ceil()) as usize;
        let pk = setup_kzg(&g1, &g2, codeword_size);

        Self::PublicParams {
            codeword_size: codeword_size,
            pk: pk,
            chunk_size: chunk_size,
        }
    }

    fn encode(data: &[u8], params: &Self::PublicParams) -> Self::EncodedData {
        let start = Instant::now();
        let rs = setup_rs2d(params.codeword_size, params.codeword_size, data.len());
        let codewords = encode_rs2d(&data, &rs);

        let result = Self::EncodedData {
            codewords: codewords,
            data_size: data.len(),
        };

        // Calculate total encoded size in bytes
        let encoded_size: usize = result.codewords.iter().map(|chunk| chunk.len()).sum();
        // Record metrics
        METRICS.with(|m| {
            let mut metrics = m.borrow_mut();
            metrics.encoding_time += start.elapsed();
            metrics.encoded_size += encoded_size;
        });

        result
    }

    fn commit(encoded: &Self::EncodedData, params: &Self::PublicParams) -> Self::Commitment {
        let start = Instant::now();

        let row_commitments = encoded
            .codewords
            .iter()
            .map(|row| {
                let poly = Polynomial {
                    coef: row.iter().map(|e| FqOrder::from_value(e.clone())).collect(),
                };
                commit_kzg(&poly, &params.pk)
            })
            .collect();

        let col_commitments = (0..encoded.codewords[0].len())
            .map(|i| {
                let col_poly = Polynomial {
                    coef: encoded
                        .codewords
                        .iter()
                        .map(|row| FqOrder::from_value(row[i].clone()))
                        .collect(),
                };
                commit_kzg(&col_poly, &params.pk)
            })
            .collect();

        let result = Self::Commitment {
            row_commitments,
            col_commitments,
        };

        let commitment_size = (result.row_commitments.len() + result.col_commitments.len())
            * std::mem::size_of::<CommitmentKZG>();

        METRICS.with(|m| {
            let mut metrics = m.borrow_mut();
            metrics.commitment_time += start.elapsed();
            metrics.commitment_size += commitment_size;
        });

        result
    }

    fn verify(
        position: &SamplePosition,
        encoded: &Self::EncodedData,
        commitment: &Self::Commitment,
        params: &Self::PublicParams,
    ) -> bool {
        let start = Instant::now();

        let poly = if position.is_row {
            Polynomial {
                coef: encoded.codewords[position.row]
                    .iter()
                    .map(|e| FqOrder::from_value(e.clone()))
                    .collect(),
            }
        } else {
            Polynomial {
                coef: encoded
                    .codewords
                    .iter()
                    .map(|r| FqOrder::from_value(r[position.col].clone()))
                    .collect(),
            }
        };

        let proof_start = Instant::now();
        // dummy value
        let x = FqOrder::from_value(5);
        let proof_kzg = open_kzg(&poly, &x, &params.pk);
        let proof_time = proof_start.elapsed();

        let sampled_commitment = if position.is_row {
            &commitment.row_commitments[position.row]
        } else {
            &commitment.col_commitments[position.col]
        };

        let result = verify_kzg(&x, sampled_commitment, &proof_kzg, &params.pk);

        // Record metrics
        METRICS.with(|m| {
            let mut metrics = m.borrow_mut();
            metrics.verification_time += start.elapsed() - proof_time; // Subtract proof time
            metrics.proof_time += proof_time;
            metrics.proof_size += std::mem::size_of::<ProofKZG>();
        });

        result
    }

    fn reconstruct(encoded: &Self::EncodedData, params: &Self::PublicParams) -> Vec<u8> {
        let _start = Instant::now();

        let rs = setup_rs2d(
            params.codeword_size,
            params.codeword_size,
            encoded.data_size,
        );
        decode_rs2d(&encoded.codewords, &rs).unwrap()
    }

    fn metrics() -> SystemMetrics {
        METRICS.with(|m| m.borrow().clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_avail_no_error() {
        // 6 * 6 to store 32 samples
        // 2x encoding
        // 12 * 12 -> 7 * 7 = 49以上欠損があるとだめ
        let params = Avail::setup(6, 2.0);

        let data: Vec<_> = (0..32).collect();
        let encoded = Avail::encode(&data, &params);
        let commit = Avail::commit(&encoded, &params);

        for i in 0..10 {
            let position0 = SamplePosition {
                row: i / 12,
                col: i % 12,
                is_row: false,
            };
            assert!(Avail::verify(&position0, &encoded, &commit, &params));
        }

        let reconstructed = Avail::reconstruct(&encoded, &params);
        for i in 0..32 {
            assert_eq!(data[i], reconstructed[i]);
        }

        METRICS.with(|metrics| {
            println!("{:#?}", *metrics.borrow());
        });

        assert!(false);

        /*
        let position1 = SamplePosition {
            row: 0,
            col: 0,
            is_row: true,
        };
        assert!(Avail::verify(&position1, &encoded, &commit, &params));

        let position2 = SamplePosition {
            row: 0,
            col: 1,
            is_row: false,
        };
        assert!(Avail::verify(&position2, &encoded, &commit, &params));

        let position3 = SamplePosition {
            row: 0,
            col: 1,
            is_row: true,
        };
        assert!(Avail::verify(&position3, &encoded, &commit, &params));

        let reconstructed = Avail::reconstruct(&encoded, &params);
        for i in 0..4 {
            assert_eq!(data[i], reconstructed[i]);
        }
        */
    }
}
