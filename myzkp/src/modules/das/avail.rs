use crate::modules::algebra::algebra::kzg::{
    commit_kzg, open_kzg, setup_kzg, verify_kzg, CommitmentKZG, ProofKZG, PublicKeyKZG,
};
use crate::modules::algebra::reedsolomon::{encode_rs2d, setup_rs2d};

pub type EncodedDataAvail = Vec<Vec<u8>>;

pub struct CommitmentAvail {
    pub row_commitments: Vec<CommitmentKZG>,
    pub col_commitments: Vec<CommitmentKZG>,
}

pub struct ProofAvail {
    pub proof_kzg: ProofKZG,
    pub is_row: bool,
}

pub struct Avail;
impl Avail {
    pub fn setup(g1: &G1Point, g2: &G2Point, width: usize, height: usize) -> PublicKeyKZG {
        setup_kzg(g1, g2, width.max(height))
    }

    pub fn encode(width: usize, height: usize, data: &[u8]) -> EncodedDataAvail {
        let rs = setup_rs2d(width, height, data.len());
        let encoded = encode_rs2d(&data, &rs);
        encoded
    }

    pub fn commit(encoded: &EncodedDataAvail, pk: &PublicKeyKZG) -> CommitmentAvail {
        let row_commits = encoded
            .iter()
            .map(|row| {
                let poly = Polynomial::from_coefficients(row);
                commit_kzg(&poly, &pk)
            })
            .collect();

        let col_commits = (0..encoded[0].len())
            .map(|i| {
                let col_poly =
                    Polynomial::from_coefficients(encoded.iter().map(|row| row[i]).collect());
                commit_kzg(&col_poly, &pk)
            })
            .collect();

        CommitmentAvail {
            row_commitments,
            col_commitments,
        }
    }

    pub fn sample(
        row_id: usize,
        col_id: usize,
        x: &FqOrder,
        encoded: &EncodedDataAvail,
        is_row: bool,
        pk: &PublicKeyKZG,
    ) -> ProofAvail {
        let poly = if is_row {
            Polynomial::from_coefficients(encoded[row_id]);
        } else {
            Polynomial::from_coefficients(encoded.iter().map(|r| r[col]).collect())
        };
        let proof_kzg = open_kzg(&poly, x, &pk);
        ProofAvail { proof_kzg, is_row }
    }

    pub fn verify(
        row_id: usize,
        col_id: usize,
        x: &FqOrder,
        commitment: &CommitmentAvail,
        proof: &ProofAvail,
        pk: &PublicKeyKZG,
    ) -> bool {
        if proof.is_row {
            verify_kzg(x, &commitment.row_commitments[row_id], proof.proof_kzg, pk)
        } else {
            verify_kzg(x, &commitment.row_commitments[col_id], proof.proof_kzg, pk)
        }
    }
}
