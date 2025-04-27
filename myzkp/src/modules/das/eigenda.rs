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
    pub expansion_factor: f64,
    pub quorums: Vec<PublicKeyKZG>,
    pub chunk_size: usize,
}

const QUORUM_COUNT: usize = 1;

pub struct EigenDA;

impl DataAvailabilitySystem for EigenDA {
    type EncodedData = EncodedDataEigenDA;
    type Commitment = CommitmentEigenDA;
    type PublicParams = PublicParamsEigenDA;

    fn setup(chunk_size: usize, expansion_factor: f64) -> Self::PublicParams {
        let g1 = BN128::generator_g1();
        let g2 = BN128::generator_g2();
        let quorums = (0..QUORUM_COUNT)
            .map(|_| setup_kzg(&g1, &g2, chunk_size)) // Standard chunk size
            .collect();

        Self::PublicParams {
            expansion_factor: expansion_factor,
            quorums: quorums,
            chunk_size: chunk_size,
        }
    }

    fn encode(data: &[u8], params: &Self::PublicParams) -> Self::EncodedData {
        let start = Instant::now();

        let codeword_size = (data.len() as f64 * params.expansion_factor.ceil()) as usize;
        let rs = setup_rs1d(codeword_size, data.len());
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

        let quorum_id: usize = 0; // dummy
        let pk = &params.quorums[quorum_id as usize];

        // commitments
        let chunk_commitments: Vec<_> = encoded
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
        let chunk_proofs: Vec<_> = encoded
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
            metrics.commitment_time += start.elapsed() - proof_time; // Subtract proof time
            metrics.proof_time += proof_time;
            metrics.commitment_size += commitment_size;
            metrics.proof_size += proof_size;
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
            metrics.verification_time += start.elapsed();
        });

        result
    }

    fn reconstruct(encoded: &Self::EncodedData, params: &Self::PublicParams) -> Vec<u8> {
        let start = Instant::now();

        let codeword_size = (encoded.data_size as f64 * params.expansion_factor.ceil()) as usize;
        let rs = setup_rs1d(codeword_size, encoded.data_size);
        let codeword = encoded.codewords.concat();
        let result = decode_rs1d(&codeword, &rs).unwrap();

        // Record metrics
        METRICS.with(|m| {
            let mut metrics = m.borrow_mut();
            metrics.reconstruction_time += start.elapsed();
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
        let params = EigenDA::setup(16, 4.0);

        let data: Vec<_> = (0..32).collect();
        let encoded = EigenDA::encode(&data, &params);
        let commit = EigenDA::commit(&encoded, &params);

        assert_eq!(encoded.codewords.len(), 8);
        assert_eq!(commit.chunk_commitments.len(), 8);
        assert_eq!(commit.chunk_proofs.len(), 8);

        for i in 0..5 {
            let position = SamplePosition {
                row: 0,
                col: i,
                is_row: false,
            };
            let isvalid = EigenDA::verify(&position, &encoded, &commit, &params);
            assert!(isvalid);
        }

        let mut received_codewords = (0..encoded.codewords.len())
            .map(|j| (0..encoded.codewords[j].len()).map(|_| 0).collect())
            .collect::<Vec<_>>();
        for i in 0..5 {
            received_codewords[i] = encoded.codewords[i].clone();
        }
        let received_encoded = EncodedDataEigenDA {
            codewords: received_codewords,
            data_size: encoded.data_size,
        };
        assert_eq!(received_encoded.codewords.len(), 8);

        let reconstructed = EigenDA::reconstruct(&received_encoded, &params);

        for i in 0..31 {
            assert_eq!(data[i], reconstructed[i]);
        }
    }
}
