use blake2::{digest::consts::U32, Blake2b, Digest};

use crate::modules::algebra::field::Field;
use crate::modules::algebra::merkle::Merkle;
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::zkstark::fiat_shamir::FiatShamirTransformer;

fn sample<F: Field>(byte_array: &[u8]) -> F {
    let mut acc: usize = 0;
    for &b in byte_array {
        acc = (acc << 8) ^ (b as usize);
    }
    F::from_value(acc)
}

fn sample_index(byte_array: &[u8], size: usize) -> usize {
    let mut acc: usize = 0;
    for &b in byte_array {
        acc = (acc << 8) ^ (b as usize);
    }
    acc % size
}

fn sample_indices(seed: &[u8], size: usize, reduced_size: usize, number: usize) -> Vec<usize> {
    assert!(
        number <= reduced_size,
        "cannot sample more indices than available in last codeword; requested: {}, available: {}",
        number,
        reduced_size
    );
    assert!(
        number <= 2 * reduced_size,
        "not enough entropy in indices wrt last codeword"
    );

    let mut indices = Vec::new();
    let mut reduced_indices = Vec::new();
    let mut counter: u64 = 0;

    while indices.len() < number {
        // Compute Blake2b hash of the concatenation of seed and counter bytes.
        let mut hasher = Blake2b::<U32>::new();
        hasher.update(seed);
        hasher.update(&counter.to_le_bytes());
        let hash_result = hasher.finalize();

        let index = sample_index(&hash_result, size);
        let reduced_index = index % reduced_size;
        counter += 1;

        // Only add the index if the reduced_index has not been seen yet.
        if !reduced_indices.contains(&reduced_index) {
            indices.push(index);
            reduced_indices.push(reduced_index);
        }
    }
    indices
}

pub struct FRI<F: Field> {
    pub offset: F,
    pub omega: F,
    pub domain_length: usize,
    pub expansion_factor: usize,
    pub num_colinearity_tests: usize,
}

pub type Codeword<F> = Vec<F>;

impl<F: Field> FRI<F> {
    fn num_rounds(&self) -> usize {
        let mut codeword_length = self.domain_length;
        let mut num_rounds = 0;
        while codeword_length > self.expansion_factor
            && 4 * self.num_colinearity_tests < codeword_length
        {
            codeword_length /= 2;
            num_rounds += 1;
        }
        num_rounds
    }

    fn eval_domain(&self) -> Vec<F> {
        (0..self.domain_length)
            .map(|i| self.offset.clone() * (self.omega.pow(i)))
            .collect()
    }

    pub fn prove(
        &self,
        codeword: &Codeword<F>,
        proof_stream: &mut FiatShamirTransformer,
    ) -> Vec<usize> {
        assert!(
            self.domain_length == codeword.len(),
            "initial codeword length does not match length of inital codeword"
        );

        // commit phase
        let codewords = self.commit(codeword, proof_stream);

        // get indices
        let top_level_indices = sample_indices(
            &proof_stream.prover_fiat_shamir(32),
            codewords[1].len(),
            codewords.last().unwrap().len(),
            self.num_colinearity_tests,
        );
        let mut indices = top_level_indices.clone();

        // query phase
        for i in 0..(codewords.len() - 1) {
            indices = indices
                .into_iter()
                .map(|idx| idx % (codewords[i].len() / 2))
                .collect::<Vec<_>>()
                .clone();
            self.query(&codewords[i], &codewords[i + 1], &indices, proof_stream);
        }

        top_level_indices
    }

    pub fn query(
        &self,
        cur_codeword: &Codeword<F>,
        next_codeword: &Codeword<F>,
        c_indices: &Vec<usize>,
        proof_stream: &mut FiatShamirTransformer,
    ) -> Vec<usize> {
        let a_indices = c_indices.clone();
        let b_indices = c_indices
            .into_iter()
            .map(|idx| idx + cur_codeword.len() / 2)
            .collect::<Vec<_>>();

        for s in 0..self.num_colinearity_tests {
            proof_stream.push(
                &bincode::serialize(&cur_codeword[a_indices[s]]).expect("Serialization failed"),
            );
            proof_stream.push(
                &bincode::serialize(&cur_codeword[b_indices[s]]).expect("Serialization failed"),
            );
            proof_stream.push(
                &bincode::serialize(&next_codeword[c_indices[s]]).expect("Serialization failed"),
            );
        }

        let serialized_cur_codeword = cur_codeword
            .into_iter()
            .map(|c| bincode::serialize(&c).expect("Serialization failed"))
            .collect::<Vec<_>>();
        let serialized_next_codeword = next_codeword
            .into_iter()
            .map(|c| bincode::serialize(&c).expect("Serialization failed"))
            .collect::<Vec<_>>();

        for s in 0..self.num_colinearity_tests {
            let path_a = Merkle::open(a_indices[s], &serialized_cur_codeword);
            let path_b = Merkle::open(b_indices[s], &serialized_cur_codeword);
            let path_c = Merkle::open(c_indices[s], &serialized_next_codeword);
            for p in path_a {
                proof_stream.push(&p);
            }
            for p in path_b {
                proof_stream.push(&p);
            }
            for p in path_c {
                proof_stream.push(&p);
            }
        }

        let mut result = a_indices.clone();
        let mut tmp = b_indices.clone();
        result.append(&mut tmp);
        result
    }

    pub fn commit(
        &self,
        initial_codeword: &Codeword<F>,
        proof_stream: &mut FiatShamirTransformer,
    ) -> Vec<Codeword<F>> {
        let one = F::one();
        let two = one.clone() + one.clone();
        let two_inverse = two.inverse();
        let mut omega = self.omega.clone();
        let mut offset = self.offset.clone();
        let mut codeword = initial_codeword.clone();
        let mut codewords = Vec::new();

        for r in 0..self.num_rounds() {
            // compute and send Merkle root
            let root = Merkle::commit(
                &codeword
                    .clone()
                    .into_iter()
                    .map(|c| bincode::serialize(&c).expect("Serialization failed"))
                    .collect::<Vec<_>>(),
            );
            proof_stream.push(&root);

            // prepare next round, if necessary
            if r == self.num_rounds() - 1 {
                break;
            }

            // get challenge
            let alpha: F = sample(&proof_stream.prover_fiat_shamir(32));

            // collect codeword
            codewords.push(codeword.clone());

            // split and fold
            codeword = (0..(codeword.len() / 2))
                .map(|i| {
                    two_inverse.clone()
                        * (one.clone() + alpha.clone() / (offset.clone() * (omega.pow(i))))
                        * codeword[i].clone()
                        + (one.clone() - alpha.clone() / (offset.clone() * (omega.pow(i))))
                            * codeword[codeword.len() / 2 + i].clone()
                })
                .collect();
            omega = omega.clone() * omega.clone();
            offset = offset.clone() * offset.clone();
        }

        // send last codeword
        let serialized_codeword = codeword
            .clone()
            .into_iter()
            .map(|c| bincode::serialize(&c).expect("Serialization failed"))
            .collect::<Vec<_>>();
        for c in serialized_codeword {
            proof_stream.push(&c);
        }

        // collect last codeword too
        codewords.push(codeword);
        codewords
    }

    /*
    pub fn verify(proof_stream: &mut FiatShamirTransformer) -> bool {
        let omega = self.omega.clone();
        let offset = self.offset.clone();

        // extract all roots and alphas
        let roots = Vec::new();
        let alphas = Vec::new();
        for r in 0..self.num_rounds() {
            roots.push(proof_stream.pull());
            alphas.push()
        }
    }*/
}
