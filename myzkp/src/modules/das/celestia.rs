use std::time::Instant;

use crate::modules::algebra::merkle::Merkle;
use crate::modules::algebra::reedsolomon::{decode_rs2d, encode_rs2d, setup_rs2d};
use crate::modules::das::utils::{DataAvailabilitySystem, SamplePosition, SystemMetrics, METRICS};

pub struct EncodedDataCelestia {
    pub codewords: Vec<Vec<Vec<u8>>>,
    pub data_size: usize, // original data size
}

pub struct CommitmentCelestia {
    pub row_roots: Vec<Vec<u8>>,
    pub col_roots: Vec<Vec<u8>>,
    pub data_root: Vec<u8>,
}

pub struct ProofCelestia {
    pub proof: Vec<Vec<u8>>,
    pub is_row: bool,
}

pub struct PublicParamsCelestia {
    pub codeword_size: usize, // the expanded side length
    pub chunk_size: usize,    // the side length when reshaping the original into a 2D array
}

pub struct Celestia;
impl DataAvailabilitySystem for Celestia {
    type EncodedData = EncodedDataCelestia;
    type Commitment = CommitmentCelestia;
    type PublicParams = PublicParamsCelestia;

    fn setup(chunk_size: usize, expansion_factor: f64, _data_size: usize) -> Self::PublicParams {
        let codeword_size = (chunk_size as f64 * expansion_factor.ceil()) as usize;

        Self::PublicParams {
            codeword_size: codeword_size,
            chunk_size: chunk_size,
        }
    }

    fn encode(data: &[u8], params: &Self::PublicParams) -> Self::EncodedData {
        let start = Instant::now();

        let rs = setup_rs2d(params.codeword_size, params.codeword_size, data.len());
        let encoded = encode_rs2d(&data, &rs);
        let reshaped: Vec<Vec<Vec<u8>>> = encoded
            .iter()
            .map(|row| {
                row.iter()
                    .map(|e| vec![e.clone()])
                    .collect::<Vec<Vec<u8>>>()
            })
            .collect();

        let result = EncodedDataCelestia {
            codewords: reshaped,
            data_size: data.len(),
        };

        let encoded_size: usize = result.codewords.iter().map(|chunk| chunk.len()).sum();
        // Record metrics
        METRICS.with(|m| {
            let mut metrics = m.borrow_mut();
            metrics.encoding_time += start.elapsed();
            metrics.encoded_size += encoded_size;
        });

        result
    }

    fn commit(encoded: &Self::EncodedData, _params: &Self::PublicParams) -> Self::Commitment {
        let start = Instant::now();

        let row_roots: Vec<Vec<u8>> = encoded
            .codewords
            .iter()
            .map(|row| Merkle::commit(&row))
            .collect();
        let col_roots: Vec<Vec<u8>> = (0..encoded.codewords[0].len())
            .map(|i| {
                let col = encoded
                    .codewords
                    .iter()
                    .map(|row| row[i].clone())
                    .collect::<Vec<Vec<u8>>>();
                Merkle::commit(&col)
            })
            .collect();

        let mut all_roots = row_roots.clone();
        all_roots.extend(col_roots.iter().cloned());
        let data_root = Merkle::commit(&all_roots);

        let commitment_size: usize = row_roots.iter().map(|v| v.len()).sum::<usize>()
            + col_roots.iter().map(|v| v.len()).sum::<usize>()
            + data_root.len();

        let result = Self::Commitment {
            row_roots,
            col_roots,
            data_root,
        };

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
        _params: &Self::PublicParams,
    ) -> bool {
        let start = Instant::now();

        let proof_start = Instant::now();
        let merkle_proof = if position.is_row {
            Merkle::open(position.col, &encoded.codewords[position.row])
        } else {
            Merkle::open(
                position.row,
                &encoded
                    .codewords
                    .iter()
                    .map(|row| row[position.col].clone())
                    .collect::<Vec<_>>(),
            )
        };
        let proof_time = proof_start.elapsed();

        let proof_size = merkle_proof.iter().map(|v| v.len()).sum::<usize>();
        let proof = ProofCelestia {
            proof: merkle_proof,
            is_row: position.is_row,
        };

        let result = if proof.is_row {
            Merkle::verify(
                &commitment.row_roots[position.row],
                position.col,
                &proof.proof,
                &encoded.codewords[position.row][position.col],
            )
        } else {
            Merkle::verify(
                &commitment.col_roots[position.col],
                position.row,
                &proof.proof,
                &encoded.codewords[position.row][position.col],
            )
        };

        // Record metrics
        METRICS.with(|m| {
            let mut metrics = m.borrow_mut();
            metrics.verification_time += start.elapsed() - proof_time; // Subtract proof time
            metrics.proof_time += proof_time;
            metrics.proof_size += proof_size;
        });

        result
    }

    fn reconstruct(encoded: &Self::EncodedData, params: &Self::PublicParams) -> Vec<u8> {
        let rs = setup_rs2d(
            params.codeword_size,
            params.codeword_size,
            encoded.data_size,
        );

        let reshaped: Vec<Vec<u8>> = encoded
            .codewords
            .iter()
            .map(|row| row.iter().map(|v| v[0]).collect::<Vec<u8>>())
            .collect();

        decode_rs2d(&reshaped, &rs).unwrap()
    }

    fn metrics() -> SystemMetrics {
        METRICS.with(|m| m.borrow().clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_celestia() {
        let params = Celestia::setup(8, 2.0, 64);

        let data: Vec<_> = (0..64).collect();
        let encoded = Celestia::encode(&data, &params);
        let commit = Celestia::commit(&encoded, &params);

        for i in 0..10 {
            let position = SamplePosition {
                row: i / 12,
                col: i % 12,
                is_row: false,
            };
            assert!(Celestia::verify(&position, &encoded, &commit, &params));
        }

        let reconstructed = Celestia::reconstruct(&encoded, &params);
        for i in 0..63 {
            assert_eq!(data[i], reconstructed[i]);
        }

        let mal_data: Vec<_> = (1..65).collect();
        let mal_encoded = Celestia::encode(&mal_data, &params);
        let mal_position = SamplePosition {
            row: 0,
            col: 0,
            is_row: false,
        };
        assert!(!Celestia::verify(
            &mal_position,
            &mal_encoded,
            &commit,
            &params
        ));
    }
}
