use std::cell::RefCell;
use std::time::Duration;

pub trait DataAvailabilitySystem {
    type EncodedData;
    type Commitment;
    type PublicParams;

    fn setup(data_size: usize, expansion_factor: f64) -> Self::PublicParams;
    fn encode(data: &[u8], params: &Self::PublicParams) -> Self::EncodedData;
    fn commit(encoded: &Self::EncodedData, params: &Self::PublicParams) -> Self::Commitment;
    fn verify(
        position: &SamplePosition,
        encoded: &Self::EncodedData,
        commitment: &Self::Commitment,
        params: &Self::PublicParams,
    ) -> bool;
    fn reconstruct(encoded: &Self::EncodedData, params: &Self::PublicParams) -> Vec<u8>;
    fn metrics() -> SystemMetrics;
}

pub struct SamplePosition {
    pub row: usize,
    pub col: usize,
    pub is_row: bool,
}

#[derive(Debug)]
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

// Thread-local metrics storage for each system
thread_local! {
    pub static METRICS: RefCell<SystemMetrics> = RefCell::new(SystemMetrics::default());
}

// Default implementation for SystemMetrics
impl Default for SystemMetrics {
    fn default() -> Self {
        SystemMetrics {
            encoding_time: Duration::default(),
            commitment_time: Duration::default(),
            proof_time: Duration::default(),
            verification_time: Duration::default(),
            reconstruction_time: Duration::default(),
            encoded_size: 0,
            commitment_size: 0,
            proof_size: 0,
        }
    }
}

// Helper function to reset metrics before benchmarking
pub fn reset_metrics() {
    METRICS.with(|m| {
        *m.borrow_mut() = SystemMetrics::default();
    });
}

// Implement Clone for SystemMetrics
impl Clone for SystemMetrics {
    fn clone(&self) -> Self {
        SystemMetrics {
            encoding_time: self.encoding_time,
            commitment_time: self.commitment_time,
            proof_time: self.proof_time,
            verification_time: self.verification_time,
            reconstruction_time: self.reconstruction_time,
            encoded_size: self.encoded_size,
            commitment_size: self.commitment_size,
            proof_size: self.proof_size,
        }
    }
}
