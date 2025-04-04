use crate::modules::algebra::curve::bn128::{FqOrder, G1Point, G2Point};
use crate::modules::algebra::kzg::{
    commit_kzg, open_kzg, setup_kzg, verify_kzg, CommitmentKZG, ProofKZG, PublicKeyKZG,
};
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::reedsolomon::{encode_rs1d, setup_rs1d};
use crate::modules::algebra::ring::Ring;

pub type EncodedDataEigenDA = Vec<Vec<u8>>;

pub struct CommitmentEigenDA {
    pub chunk_commitments: Vec<CommitmentKZG>, // Per-chunk commitments
    pub quorum_id: u32, // Multiple quorums supported
}

pub struct ProofEigenDA {
    pub proof_kzg: ProofKZG,
    pub chunk_index: usize,
    pub quorum_id: u32,
}

pub struct EigenDisperser {
    pub expansion_factor: usize,
    pub quorums: Vec<PublicKeyKZG>,
}

impl EigenDisperser {
    pub fn new(g1: &G1Point, g2: &G2Point, chunk_size: usize, expansion_factor: usize, quorum_count: usize) -> Self {
        let quorums = (0..quorum_count)
            .map(|_| setup_kzg(g1, g2, chunk_size)) // Standard chunk size
            .collect();
            
        EigenDisperser {
            expansion_factor,
            quorums,
        }
    }

    pub fn encode(&self, data: &[u8], chunk_size: usize) -> EncodedDataEigenDA {
        let rs = setup_rs1d(self.expansion_factor, data.len());
        let encoded = encode_rs1d(data, &rs);
        encoded.chunks(chunk_size) // Fixed chunk size
            .map(|c| c.to_vec())
            .collect()
    }

    pub fn commit(&self, encoded: &EncodedDataEigenDA, quorum_id: u32) -> CommitmentEigenDA {
        let pk = &self.quorums[quorum_id as usize];
        let chunk_commitments = encoded
            .iter()
            .map(|chunk| {
                let poly = Polynomial {
                    coef: chunk.iter().map(|e| FqOrder::from_value(e.clone())).collect(),
                };
                commit_kzg(&poly, pk)
            })
            .collect();

        CommitmentEigenDA {
            chunk_commitments,
            quorum_id,
        }
    }
}

pub struct EigenDA;
impl EigenDA {
    pub fn sample(
        chunk_index: usize,
        x: &FqOrder,
        encoded: &EncodedDataEigenDA,
        quorum_id: u32,
        disperser: &EigenDisperser,
    ) -> ProofEigenDA {
        let pk = &disperser.quorums[quorum_id as usize];
        let poly = Polynomial {
            coef: encoded[chunk_index]
                .iter()
                .map(|e| FqOrder::from_value(e.clone()))
                .collect(),
        };
        let proof_kzg = open_kzg(&poly, x, pk);
        
        ProofEigenDA {
            proof_kzg,
            chunk_index,
            quorum_id,
        }
    }

    pub fn verify(
        x: &FqOrder,
        commitment: &CommitmentEigenDA,
        proof: &ProofEigenDA,
        disperser: &EigenDisperser,
    ) -> bool {
        let pk = &disperser.quorums[proof.quorum_id as usize];
        verify_kzg(
            x,
            &commitment.chunk_commitments[proof.chunk_index],
            &proof.proof_kzg,
            pk,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modules::algebra::curve::bn128::BN128;

    #[test]
    fn test_eigenda_flow() {
        let g1 = BN128::generator_g1();
        let g2 = BN128::generator_g2();
        
        // Initialize disperser with 4.5x redundancy and 3 quorums
        let disperser = EigenDisperser::new(&g1, &g2, 2, 9, 3); 
        
        let data = vec![1, 2, 3, 4];
        let encoded = disperser.encode(&data, 2);
        // assert_eq!(encoded.len(), 9); // 4.5x expansion
        
        let commitment = disperser.commit(&encoded, 0); // Commit to quorum 0
        
        let z = FqOrder::from_value(5);
        let proof = EigenDA::sample(3, &z, &encoded, 0, &disperser);
        assert!(EigenDA::verify(&z, &commitment, &proof, &disperser));
    }
}