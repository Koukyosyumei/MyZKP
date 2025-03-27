use crate::modules::algebra::curve::bn128::{FqOrder, G1Point, G2Point};
use crate::modules::algebra::kzg::{
    commit_kzg, open_kzg, setup_kzg, verify_kzg, CommitmentKZG, ProofKZG, PublicKeyKZG,
};
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::reedsolomon::{encode_rs2d, setup_rs2d};
use crate::modules::algebra::ring::Ring;

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
        let row_commitments = encoded
            .iter()
            .map(|row| {
                let poly = Polynomial {
                    coef: row.iter().map(|e| FqOrder::from_value(e.clone())).collect(),
                };
                commit_kzg(&poly, &pk)
            })
            .collect();

        let col_commitments = (0..encoded[0].len())
            .map(|i| {
                let col_poly = Polynomial {
                    coef: encoded
                        .iter()
                        .map(|row| FqOrder::from_value(row[i].clone()))
                        .collect(),
                };
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
            Polynomial {
                coef: encoded[row_id]
                    .iter()
                    .map(|e| FqOrder::from_value(e.clone()))
                    .collect(),
            }
        } else {
            Polynomial {
                coef: encoded
                    .iter()
                    .map(|r| FqOrder::from_value(r[col_id].clone()))
                    .collect(),
            }
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
            verify_kzg(x, &commitment.row_commitments[row_id], &proof.proof_kzg, pk)
        } else {
            verify_kzg(x, &commitment.row_commitments[col_id], &proof.proof_kzg, pk)
        }
    }
}
