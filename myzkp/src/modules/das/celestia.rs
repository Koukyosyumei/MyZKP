use crate::modules::algebra::merkle::Merkle;
use crate::modules::algebra::reedsolomon::{decode_rs2d, encode_rs2d, setup_rs2d};

pub type EncodedDataCelestia = Vec<Vec<Vec<u8>>>;

pub struct CommitmentCelestia {
    pub row_roots: Vec<Vec<u8>>,
    pub col_roots: Vec<Vec<u8>>,
}

pub struct ProofCelestia {
    pub row_proof: Vec<Vec<u8>>,
    pub col_proof: Vec<Vec<u8>>,
}

pub struct Celestia;
impl Celestia {
    pub fn encode(
        width: usize,
        height: usize,
        data: &[u8],
    ) -> (EncodedDataCelestia, CommitmentCelestia) {
        let rs = setup_rs2d(width, height, data.len());
        let encoded = encode_rs2d(&data, &rs);
        let reshaped: EncodedDataCelestia = encoded
            .iter()
            .map(|row| {
                row.iter()
                    .map(|e| vec![e.clone()])
                    .collect::<Vec<Vec<u8>>>()
            })
            .collect();

        let row_roots: Vec<Vec<u8>> = reshaped.iter().map(|row| Merkle::commit(&row)).collect();
        let col_roots: Vec<Vec<u8>> = (0..reshaped[0].len())
            .map(|i| {
                let col = reshaped
                    .iter()
                    .map(|row| row[i].clone())
                    .collect::<Vec<Vec<u8>>>();
                Merkle::commit(&col)
            })
            .collect();

        (
            reshaped,
            CommitmentCelestia {
                row_roots,
                col_roots,
            },
        )
    }

    pub fn sample(row_id: usize, col_id: usize, encoded: &EncodedDataCelestia) -> ProofCelestia {
        let row_proof = Merkle::open(col_id, &encoded[row_id]);
        let col_proof = Merkle::open(
            row_id,
            &encoded
                .iter()
                .map(|row| row[col_id].clone())
                .collect::<Vec<_>>(),
        );

        ProofCelestia {
            row_proof,
            col_proof,
        }
    }

    pub fn verify(
        row_id: usize,
        col_id: usize,
        chunk: &[u8],
        commitment: &CommitmentCelestia,
        proof: &ProofCelestia,
    ) -> bool {
        Merkle::verify(
            &commitment.row_roots[row_id],
            col_id,
            &proof.row_proof,
            chunk,
        ) && Merkle::verify(
            &commitment.row_roots[col_id],
            row_id,
            &proof.col_proof,
            chunk,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_celestia_no_error() {
        let data = vec![1, 2, 3, 4];
        let (encoded, commitment) = Celestia::encode(2, 2, &data);
        let proof = Celestia::sample(0, 0, &encoded);
        assert!(Celestia::verify(0, 0, &encoded[0][0], &commitment, &proof));
    }
}
