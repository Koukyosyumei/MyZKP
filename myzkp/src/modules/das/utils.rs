pub trait DataAvailabilitySystem {
    type EncodedData;
    type Commitment;
    type PublicParams;

    fn setup(data_size: usize, expansion_factor: f64) -> Self::PublicParams;
    fn encode(data: &[u8], params: &Self::PublicParams) -> Self::EncodedData;
    fn commit(encoded: &Self::EncodedData, params: &Self::PublicParams) -> Self::Commitment;
    fn verify(
        position: &SamplePosition,
        commitment: &Self::Commitment,
        params: &Self::PublicParams,
    ) -> bool;
    //fn reconstruct(encoded: &Self::EncodedData, params: &Self::PublicParams) -> Vec<u8>;
    //fn metrics() -> SystemMetrics;
}

pub struct SamplePosition {
    pub row: usize,
    pub col: usize,
}

pub struct SystemMetrics {
    pub encoding_time: std::time::Duration,
    pub commitment_time: std::time::Duration,
    pub proof_time: std::time::Duration,
    pub verification_time: std::time::Duration,
    pub reconstruction_time: std::time::Duration,
    pub encoded_size: usize,
    pub commitment_size: usize,
    pub proof_size: usize,
}
