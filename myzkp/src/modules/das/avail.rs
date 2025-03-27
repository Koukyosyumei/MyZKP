use crate::modules::algebra::algebra::kzg::{
    commit_kzg, open_kzg, setup_kzg, verify_kzg, CommitmentKZG, ProofKZG, PublicKeyKZG,
};
use crate::modules::algebra::reedsolomon::{encode_rs2d, setup_rs2d};

pub struct Avail;

impl Avail {
    pub fn setup(g1: &G1Point, g2: &G2Point, width: usize, height: usize) -> PublicKeyKZG {
        setup_kzg(g1, g2, width.max(height))
    }
}
