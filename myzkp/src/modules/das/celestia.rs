use crate::modules::algebra::merkle::Merkle;
use crate::modules::algebra::reedsolomon::{decode_rs2d, encode_rs2d, setup_rs2d};
use crate::modules::das::utils::{DataAvailabilitySystem, SamplePosition};

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

    fn setup(chunk_size: usize, expansion_factor: f64) -> Self::PublicParams {
        let codeword_size = (chunk_size as f64 * expansion_factor.ceil()) as usize;

        Self::PublicParams {
            codeword_size: codeword_size,
            chunk_size: chunk_size,
        }
    }

    fn encode(data: &[u8], params: &Self::PublicParams) -> Self::EncodedData {
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

        EncodedDataCelestia {
            codewords: reshaped,
            data_size: data.len(),
        }
    }

    fn commit(encoded: &Self::EncodedData, params: &Self::PublicParams) -> Self::Commitment {
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

        Self::Commitment {
            row_roots,
            col_roots,
            data_root,
        }
    }

    fn verify(
        position: &SamplePosition,
        encoded: &Self::EncodedData,
        commitment: &Self::Commitment,
        params: &Self::PublicParams,
    ) -> bool {
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
        let proof = ProofCelestia {
            proof: merkle_proof,
            is_row: position.is_row,
        };

        if proof.is_row {
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
        }
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_celestia_no_error() {
        let params = Celestia::setup(2, 2.0);

        let data = vec![1, 2, 3, 4];
        let encoded = Celestia::encode(&data, &params);
        let commit = Celestia::commit(&encoded, &params);

        let position0 = SamplePosition {
            row: 0,
            col: 0,
            is_row: false,
        };
        assert!(Celestia::verify(&position0, &encoded, &commit, &params));

        let position1 = SamplePosition {
            row: 0,
            col: 0,
            is_row: true,
        };
        assert!(Celestia::verify(&position1, &encoded, &commit, &params));

        let position2 = SamplePosition {
            row: 0,
            col: 1,
            is_row: false,
        };
        assert!(Celestia::verify(&position2, &encoded, &commit, &params));

        let position3 = SamplePosition {
            row: 0,
            col: 1,
            is_row: true,
        };
        assert!(Celestia::verify(&position3, &encoded, &commit, &params));

        let reconstructed = Celestia::reconstruct(&encoded, &params);
        for i in 0..4 {
            assert_eq!(data[i], reconstructed[i]);
        }
    }
}
