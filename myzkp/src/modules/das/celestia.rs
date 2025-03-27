use crate::modules::algebra::merkle::Merkle;
use crate::modules::algebra::reedsolomon::{decode_rs2d, encode_rs2d, setup_rs2d};

pub type EncodedDataCelestia = Vec<Vec<Vec<u8>>>;

pub struct CommitmentCelestia {
    pub row_roots: Vec<Vec<u8>>,
    pub col_roots: Vec<Vec<u8>>,
    pub data_root: Vec<u8>,
}

pub struct ProofCelestia {
    pub proof: Vec<Vec<u8>>,
    pub is_left: bool,
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

        let mut all_roots = row_roots.clone();
        all_roots.extend(col_roots.iter().cloned());
        let data_root = Merkle::commit(&all_roots);

        (
            reshaped,
            CommitmentCelestia {
                row_roots,
                col_roots,
                data_root,
            },
        )
    }

    pub fn sample(
        row_id: usize,
        col_id: usize,
        encoded: &EncodedDataCelestia,
        is_left: bool,
    ) -> ProofCelestia {
        let proof = if is_left {
            Merkle::open(col_id, &encoded[row_id])
        } else {
            Merkle::open(
                row_id,
                &encoded
                    .iter()
                    .map(|row| row[col_id].clone())
                    .collect::<Vec<_>>(),
            )
        };

        ProofCelestia { proof, is_left }
    }

    pub fn verify(
        row_id: usize,
        col_id: usize,
        chunk: &[u8],
        commitment: &CommitmentCelestia,
        proof: &ProofCelestia,
    ) -> bool {
        if proof.is_left {
            Merkle::verify(&commitment.row_roots[row_id], col_id, &proof.proof, chunk)
        } else {
            Merkle::verify(&commitment.col_roots[col_id], row_id, &proof.proof, chunk)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_celestia_no_error() {
        let data = vec![1, 2, 3, 4];
        let (encoded, commitment) = Celestia::encode(4, 4, &data);
        let proof_00_left = Celestia::sample(0, 0, &encoded, true);
        assert!(Celestia::verify(
            0,
            0,
            &encoded[0][0],
            &commitment,
            &proof_00_left
        ));
        let proof_00_right = Celestia::sample(0, 0, &encoded, false);
        assert!(Celestia::verify(
            0,
            0,
            &encoded[0][0],
            &commitment,
            &proof_00_right
        ));
        let proof_01_left = Celestia::sample(0, 1, &encoded, true);
        assert!(Celestia::verify(
            0,
            1,
            &encoded[0][1],
            &commitment,
            &proof_01_left
        ));
        let proof_01_right = Celestia::sample(0, 1, &encoded, false);
        assert!(Celestia::verify(
            0,
            1,
            &encoded[0][1],
            &commitment,
            &proof_01_right
        ));
        assert!(!Celestia::verify(
            0,
            1,
            &encoded[1][0],
            &commitment,
            &proof_01_left
        ));
        assert!(!Celestia::verify(
            0,
            1,
            &encoded[1][0],
            &commitment,
            &proof_01_right
        ));
    }
}
