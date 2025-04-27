use std::time::Instant;

use crate::modules::algebra::curve::bn128::FqOrder;
use crate::modules::algebra::curve::bn128::BN128;
use crate::modules::algebra::kzg::{
    commit_kzg, open_kzg, setup_kzg, verify_kzg, CommitmentKZG, ProofKZG, PublicKeyKZG,
};
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::reedsolomon::{decode_rs1d, encode_rs1d, setup_rs1d};
use crate::modules::algebra::ring::Ring;
use crate::modules::das::utils::{DataAvailabilitySystem, SamplePosition, SystemMetrics, METRICS};

pub struct EncodedDataAvail {
    pub codewords: Vec<Vec<u8>>,
    pub data_size: usize, // original data size
}

pub struct CommitmentAvail {
    pub commitments: Vec<CommitmentKZG>,
}

pub struct PublicParamsAvail {
    pub expansion_factor: f64,
    pub pk: PublicKeyKZG,
    pub chunk_size: usize,
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
            expansion_factor: expansion_factor,
            pk: pk,
            chunk_size: chunk_size,
        }
    }

    fn encode(data: &[u8], params: &Self::PublicParams) -> Self::EncodedData {
        let start = Instant::now();

        let codeword_size = (params.chunk_size as f64 * params.expansion_factor.ceil()) as usize;
        let rs = setup_rs1d(codeword_size, params.chunk_size);

        let codewords: Vec<Vec<u8>> = data
            .chunks(params.chunk_size)
            .map(|c| {
                let mut padded = c.to_vec().clone();
                padded.resize(params.chunk_size, 0);
                encode_rs1d(&padded, &rs)
            })
            .collect();

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

        let commitments = encoded
            .codewords
            .iter()
            .map(|c| {
                let poly = Polynomial {
                    coef: c.iter().map(|e| FqOrder::from_value(e.clone())).collect(),
                };
                commit_kzg(&poly, &params.pk)
            })
            .collect();

        let result = Self::Commitment { commitments };

        let commitment_size = result.commitments.len() * std::mem::size_of::<CommitmentKZG>();

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

        let poly = Polynomial {
            coef: encoded.codewords[position.col]
                .iter()
                .map(|r| FqOrder::from_value(r.clone()))
                .collect(),
        };

        let proof_start = Instant::now();
        // dummy value
        let x = FqOrder::from_value(5);
        let proof_kzg = open_kzg(&poly, &x, &params.pk);
        let proof_time = proof_start.elapsed();

        let sampled_commitment = &commitment.commitments[position.col];

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

        let codeword_size = (params.chunk_size as f64 * params.expansion_factor.ceil()) as usize;
        let rs = setup_rs1d(codeword_size, params.chunk_size);

        let decodeds: Vec<_> = encoded
            .codewords
            .iter()
            .map(|c| decode_rs1d(&c, &rs).unwrap())
            .collect();
        let mut concatenated = decodeds.concat();
        concatenated.resize(encoded.data_size, 0);
        concatenated
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
        let params = Avail::setup(8, 2.0);

        let data: Vec<_> = (0..32).collect();
        let encoded = Avail::encode(&data, &params);
        let commit = Avail::commit(&encoded, &params);

        for i in 0..3 {
            let position0 = SamplePosition {
                row: 0,
                col: i % 12,
                is_row: false,
            };
            assert!(Avail::verify(&position0, &encoded, &commit, &params));
        }

        let reconstructed = Avail::reconstruct(&encoded, &params);
        for i in 0..32 {
            assert_eq!(data[i], reconstructed[i]);
        }
    }
}
