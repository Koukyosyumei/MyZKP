use std::time::Instant;

use crate::modules::algebra::curve::bn128::FqOrder;
use crate::modules::algebra::curve::bn128::BN128;
use crate::modules::algebra::kzg::{
    commit_kzg, open_kzg, setup_kzg, verify_kzg, CommitmentKZG, ProofKZG, PublicKeyKZG,
};
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::reedsolomon::{decode_rs1d, encode_rs1d, setup_rs1d};
use crate::modules::algebra::ring::Ring;
use crate::modules::das::utils::SamplePosition;
use crate::modules::das::utils::{DataAvailabilitySystem, SystemMetrics, METRICS};

pub struct EncodedDataEigenDA {
    pub codewords: Vec<Vec<u8>>,
    pub data_size: usize, // original data size
}

pub struct CommitmentEigenDA {
    pub chunk_commitments: Vec<CommitmentKZG>, // Per-chunk commitments
    pub chunk_proofs: Vec<ProofKZG>,
    pub quorum_id: u32, // Multiple quorums supported
}

pub struct PublicParamsEigenDA {
    pub codeword_size: usize,
    pub quorums: Vec<PublicKeyKZG>,
    pub chunk_size: usize,
}

const QUORUM_COUNT: usize = 2;

pub struct EigenDA;

impl DataAvailabilitySystem for EigenDA {
    type EncodedData = EncodedDataEigenDA;
    type Commitment = CommitmentEigenDA;
    type PublicParams = PublicParamsEigenDA;

    fn setup(chunk_size: usize, expansion_factor: f64) -> Self::PublicParams {
        let g1 = BN128::generator_g1();
        let g2 = BN128::generator_g2();
        let codeword_size = (chunk_size as f64 * expansion_factor.ceil()) as usize;
        let quorums = (0..QUORUM_COUNT)
            .map(|_| setup_kzg(&g1, &g2, chunk_size)) // Standard chunk size
            .collect();

        Self::PublicParams {
            codeword_size: codeword_size,
            quorums: quorums,
            chunk_size: chunk_size,
        }
    }

    fn encode(data: &[u8], params: &Self::PublicParams) -> Self::EncodedData {
        let start = Instant::now();

        let rs = setup_rs1d(params.codeword_size, data.len());
        let encoded = encode_rs1d(data, &rs);
        let codewords: Vec<Vec<u8>> = encoded
            .chunks(params.chunk_size) // Fixed chunk size
            .map(|c| c.to_vec())
            .collect();

        let result = Self::EncodedData {
            codewords: codewords,
            data_size: data.len(),
        };

        // Calculate total encoded size in bytes
        let encoded_size = result.codewords.iter().map(|chunk| chunk.len()).sum();
        // Record metrics
        METRICS.with(|m| {
            let mut metrics = m.borrow_mut();
            metrics.encoding_time = start.elapsed();
            metrics.encoded_size = encoded_size;
        });

        result
    }

    fn commit(encoded: &Self::EncodedData, params: &Self::PublicParams) -> Self::Commitment {
        let start = Instant::now();

        let quorum_id: usize = 0; // dummy
        let pk = &params.quorums[quorum_id as usize];

        // commitments
        let chunk_commitments = encoded
            .codewords
            .iter()
            .map(|chunk| {
                let poly = Polynomial {
                    coef: chunk
                        .iter()
                        .map(|e| FqOrder::from_value(e.clone()))
                        .collect(),
                };
                commit_kzg(&poly, pk)
            })
            .collect();

        // Start proof time measurement
        let proof_start = Instant::now();

        // dummy value
        let x = FqOrder::from_value(5);
        // proofs
        let chunk_proofs = encoded
            .codewords
            .iter()
            .map(|chunk| {
                let poly = Polynomial {
                    coef: chunk
                        .iter()
                        .map(|e| FqOrder::from_value(e.clone()))
                        .collect(),
                };
                open_kzg(&poly, &x, pk)
            })
            .collect();
        let proof_time = proof_start.elapsed();

        let result = Self::Commitment {
            chunk_commitments,
            chunk_proofs,
            quorum_id: quorum_id.try_into().unwrap(),
        };

        // Calculate sizes (using std::mem::size_of_val would be more accurate)
        let commitment_size = result.chunk_commitments.len() * std::mem::size_of::<CommitmentKZG>();
        let proof_size = result.chunk_proofs.len() * std::mem::size_of::<ProofKZG>();

        // Record metrics
        METRICS.with(|m| {
            let mut metrics = m.borrow_mut();
            metrics.commitment_time = start.elapsed() - proof_time; // Subtract proof time
            metrics.proof_time = proof_time;
            metrics.commitment_size = commitment_size;
            metrics.proof_size = proof_size;
        });

        result
    }

    fn verify(
        position: &SamplePosition,
        _encoded: &Self::EncodedData,
        commitment: &Self::Commitment,
        params: &Self::PublicParams,
    ) -> bool {
        let start = Instant::now();

        let pk = &params.quorums[0];
        // dummy value
        let x = FqOrder::from_value(5);

        let result = verify_kzg(
            &x,
            &commitment.chunk_commitments[position.col],
            &commitment.chunk_proofs[position.col],
            pk,
        );

        // Record metrics
        METRICS.with(|m| {
            let mut metrics = m.borrow_mut();
            metrics.verification_time = start.elapsed();
        });

        result
    }

    fn reconstruct(encoded: &Self::EncodedData, params: &Self::PublicParams) -> Vec<u8> {
        let start = Instant::now();

        let rs = setup_rs1d(params.codeword_size, encoded.data_size);
        let codeword = encoded.codewords.concat();
        let result = decode_rs1d(&codeword, &rs).unwrap();

        // Record metrics
        METRICS.with(|m| {
            let mut metrics = m.borrow_mut();
            metrics.reconstruction_time = start.elapsed();
        });

        result
    }

    fn metrics() -> SystemMetrics {
        METRICS.with(|m| m.borrow().clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eigenda_flow() {
        let params = EigenDA::setup(2, 4.5);

        let data = vec![1, 2, 3, 4];
        let encoded = EigenDA::encode(&data, &params);
        let commit = EigenDA::commit(&encoded, &params);

        let position = SamplePosition {
            row: 0,
            col: 3,
            is_row: false,
        };
        let isvalid = EigenDA::verify(&position, &encoded, &commit, &params);
        assert!(isvalid);

        let reconstructed = EigenDA::reconstruct(&encoded, &params);
        for i in 0..4 {
            assert_eq!(data[i], reconstructed[i]);
        }
    }
}
