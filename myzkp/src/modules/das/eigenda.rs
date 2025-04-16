use crate::modules::algebra::curve::bn128::BN128;
use crate::modules::algebra::curve::bn128::{FqOrder, G1Point, G2Point};
use crate::modules::algebra::kzg::{
    commit_kzg, open_kzg, setup_kzg, verify_kzg, CommitmentKZG, ProofKZG, PublicKeyKZG,
};
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::reedsolomon::{encode_rs1d, setup_rs1d};
use crate::modules::algebra::ring::Ring;
use crate::modules::das::utils::DataAvailabilitySystem;
use crate::modules::das::utils::SamplePosition;

pub type EncodedDataEigenDA = Vec<Vec<u8>>;

pub struct CommitmentEigenDA {
    pub chunk_commitments: Vec<CommitmentKZG>, // Per-chunk commitments
    pub chunk_proofs: Vec<ProofKZG>,
    pub quorum_id: u32, // Multiple quorums supported
}

pub struct PublicParamsEigenDA {
    pub expansion_factor: usize,
    pub quorums: Vec<PublicKeyKZG>,
    pub chunk_size: usize,
}

struct EigenDA;

impl DataAvailabilitySystem for EigenDA {
    type EncodedData = EncodedDataEigenDA;
    type Commitment = CommitmentEigenDA;
    type PublicParams = PublicParamsEigenDA;

    fn setup(_data_size: usize, expansion_factor: f64) -> Self::PublicParams {
        let g1 = BN128::generator_g1();
        let g2 = BN128::generator_g2();

        // Calculate appropriate chunk size based on data size
        // let optimal_chunk_size = (data_size as f64).sqrt().ceil() as usize;
        let chunk_size = 2; // optimal_chunk_size.max(16).min(1024); // Reasonable bounds

        // Convert expansion factor to integer for RS encoding
        let expansion_factor_int = (chunk_size as f64 * expansion_factor.ceil()) as usize;

        let quorum_count = 3; // Default, could be made configurable
        let quorums = (0..quorum_count)
            .map(|_| setup_kzg(&g1, &g2, chunk_size)) // Standard chunk size
            .collect();

        PublicParamsEigenDA {
            expansion_factor: expansion_factor_int,
            quorums: quorums,
            chunk_size: chunk_size,
        }
    }

    fn encode(data: &[u8], params: &Self::PublicParams) -> Self::EncodedData {
        let rs = setup_rs1d(params.expansion_factor, data.len());
        let encoded = encode_rs1d(data, &rs);
        encoded
            .chunks(params.chunk_size) // Fixed chunk size
            .map(|c| c.to_vec())
            .collect()
    }

    fn commit(encoded: &Self::EncodedData, params: &Self::PublicParams) -> Self::Commitment {
        let quorum_id: usize = 0; // dummy
        let pk = &params.quorums[quorum_id as usize];

        // commitments
        let chunk_commitments = encoded
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

        // dummy value
        let x = FqOrder::from_value(5);

        // proofs
        let chunk_proofs = encoded
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

        Self::Commitment {
            chunk_commitments,
            chunk_proofs,
            quorum_id: quorum_id.try_into().unwrap(),
        }
    }

    fn verify(
        position: &SamplePosition,
        commitment: &Self::Commitment,
        params: &Self::PublicParams,
    ) -> bool {
        let pk = &params.quorums[0];

        // dummy value
        let x = FqOrder::from_value(5);

        verify_kzg(
            &x,
            &commitment.chunk_commitments[position.col],
            &commitment.chunk_proofs[position.col],
            pk,
        )
    }

    //fn reconstruct(encoded: &Self::EncodedData, params: &Self::PublicParams) -> Vec<u8> {
    //    todo!()
    //}

    /*
        type EncodedData;
        type Commitment;
        type Proof;
        type PublicParams;

        fn setup(data_size: usize, expansion_factor: f64) -> Self::PublicParams;
        fn encode(data: &[u8], params: &Self::PublicParams) -> Self::EncodedData;
        fn commit(encoded: &Self::EncodedData, params: &Self::PublicParams) -> Self::Commitment;
        fn sample(
            position: &SamplePosition,
            encoded: &Self::EncodedData,
            params: &Self::PublicParams,
        ) -> Self::Proof;
        fn verify(
            position: &SamplePosition,
            commitment: &Self::Commitment,
            proof: &Self::Proof,
            params: &Self::PublicParams,
        ) -> bool;
        fn reconstruct(encoded: &Self::EncodedData, params: &Self::PublicParams) -> Vec<u8>;
        fn metrics() -> SystemMetrics;
    */
}

/*
impl EigenDisperser {
    pub fn new(
        g1: &G1Point,
        g2: &G2Point,
        chunk_size: usize,
        expansion_factor: usize,
        quorum_count: usize,
    ) -> Self {
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
        encoded
            .chunks(chunk_size) // Fixed chunk size
            .map(|c| c.to_vec())
            .collect()
    }

    pub fn commit(&self, encoded: &EncodedDataEigenDA, quorum_id: u32) -> CommitmentEigenDA {
        let pk = &self.quorums[quorum_id as usize];
        let chunk_commitments = encoded
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

        CommitmentEigenDA {
            chunk_commitments,
            quorum_id,
        }
    }

    pub fn prove(
        &self,
        encoded: &EncodedDataEigenDA,
        quorum_id: u32,
        chunk_index: usize,
        x: &FqOrder,
    ) -> ProofEigenDA {
        let pk = &self.quorums[quorum_id as usize];
        let poly = Polynomial {
            coef: encoded[chunk_index]
                .iter()
                .map(|e| FqOrder::from_value(*e))
                .collect(),
        };
        let proof_kzg = open_kzg(&poly, x, pk);
        ProofEigenDA {
            proof_kzg,
            chunk_index,
            quorum_id,
        }
    }
}

pub struct EigenDA;
impl EigenDA {
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

    */

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modules::algebra::curve::bn128::BN128;

    #[test]
    fn test_eigenda_flow() {
        let params = EigenDA::setup(4, 4.5);

        let data = vec![1, 2, 3, 4];
        let encode = EigenDA::encode(&data, &params);
        let commit = EigenDA::commit(&encode, &params);

        let position = SamplePosition { row: 0, col: 3 };
        let isvalid = EigenDA::verify(&position, &commit, &params);
        //let g1 = BN128::generator_g1();
        //let g2 = BN128::generator_g2();

        /*
        // Initialize disperser with 4.5x redundancy and 3 quorums
        let disperser = EigenDisperser::new(&g1, &g2, 2, 9, 3);

        let data = vec![1, 2, 3, 4];
        let encoded = disperser.encode(&data, 2);
        // assert_eq!(encoded.len(), 9); // 4.5x expansion

        let commitment = disperser.commit(&encoded, 0); // Commit to quorum 0
        let z = FqOrder::from_value(5);
        let proof = disperser.prove(&encoded, 0, 3, &z);
        assert!(EigenDA::verify(&z, &commitment, &proof, &disperser));
        */
    }
}
