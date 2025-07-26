use std::ops::{Shl, Shr};
use std::str::FromStr;

use blake2::{digest::consts::U32, Blake2b, Digest};
use lazy_static::lazy_static;
use num_bigint::BigInt;
use num_traits::{One, Zero};
use paste::paste;
use serde::{Deserialize, Serialize};

use crate::modules::algebra::efield::{ExtendedFieldElement, IrreduciblePoly};
use crate::modules::algebra::field::{Field, FiniteFieldElement, ModulusValue};
use crate::modules::algebra::merkle::{Merkle, MerklePath, MerkleRoot};
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::ring::Ring;
use crate::modules::zkstark::fiat_shamir::FiatShamirTransformer;
use crate::{define_extension_field, define_myzkp_modulus_type};

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

pub struct FriProof<F: Field> {
    pub top_level_indices: Vec<usize>,
    pub last_codeword: Vec<F>,
    pub merkle_roots: Vec<MerkleRoot>,
    pub revealed_layers: Vec<FriQueryLayer<F>>,
}

pub struct FriQueryLayer<F: Field> {
    pub a: (Vec<F>, Vec<MerklePath>),
    pub b: (Vec<F>, Vec<MerklePath>),
    pub c: (Vec<F>, Vec<MerklePath>),
}

impl<F: Field> FRI<F> {
    pub fn num_rounds(&self) -> usize {
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

    pub fn eval_domain(&self) -> Vec<F> {
        (0..self.domain_length)
            .map(|i| self.offset.clone() * (self.omega.pow(i)))
            .collect()
    }

    pub fn prove(&self, initial_codeword: &Codeword<F>) -> FriProof<F> {
        let mut proof_stream = FiatShamirTransformer::new();

        assert!(
            self.domain_length == initial_codeword.len(),
            "initial codeword length does not match length of inital codeword"
        );

        let (codewords, roots) = self.commit(initial_codeword, &mut proof_stream);

        // get indices
        let top_level_indices = sample_indices(
            &proof_stream.prover_fiat_shamir(32),
            codewords[1].len(),
            codewords.last().unwrap().len(),
            self.num_colinearity_tests,
        );
        let mut indices = top_level_indices.clone();

        // query phase
        let mut revealed_items = Vec::new();
        for i in 0..(codewords.len() - 1) {
            indices = indices
                .into_iter()
                .map(|idx| idx % (codewords[i].len() / 2))
                .collect::<Vec<_>>()
                .clone();
            let tmp = self.reveal(&codewords[i], &codewords[i + 1], &indices);
            revealed_items.push(tmp);
        }

        FriProof {
            top_level_indices: top_level_indices,
            last_codeword: codewords.last().unwrap().to_vec(),
            merkle_roots: roots,
            revealed_layers: revealed_items,
        }
    }

    pub fn commit(
        &self,
        initial_codeword: &Codeword<F>,
        proof_stream: &mut FiatShamirTransformer,
    ) -> (Vec<Codeword<F>>, Vec<MerkleRoot>) {
        let one = F::one();
        let two = one.clone() + one.clone();
        let two_inverse = two.inverse();
        let mut omega = self.omega.clone();
        let mut offset = self.offset.clone();
        let mut codeword = initial_codeword.clone();
        let mut codewords = Vec::new();

        let mut roots = Vec::new();
        for r in 0..self.num_rounds() {
            // compute and send Merkle root
            let root = Merkle::commit(
                &codeword
                    .clone()
                    .into_iter()
                    .map(|c| bincode::serialize(&c).expect("Serialization failed"))
                    .collect::<Vec<_>>(),
            );
            roots.push(root.clone());
            proof_stream.push(&vec![root]);

            // prepare next round, if necessary
            if r == self.num_rounds() - 1 {
                break;
            }

            // get challenge
            let alpha = F::sample(&proof_stream.prover_fiat_shamir(32));

            // collect codeword
            codewords.push(codeword.clone());

            // split and fold
            codeword = (0..(codeword.len() / 2))
                .map(|i| {
                    two_inverse
                        .mul_ref(
                            &((one.add_ref(&(alpha.div_ref(&offset.mul_ref(&omega.pow(i))))))
                                .mul_ref(&codeword[i])
                                + (one.sub_ref(&(alpha.div_ref(&offset.mul_ref(&omega.pow(i))))))
                                    .mul_ref(&codeword[codeword.len() / 2 + i])),
                        )
                        .sanitize()
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
        proof_stream.push(&serialized_codeword);
        codewords.push(codeword);

        (codewords, roots)
    }

    pub fn reveal(
        &self,
        cur_codeword: &Codeword<F>,
        next_codeword: &Codeword<F>,
        c_indices: &Vec<usize>,
    ) -> FriQueryLayer<F> {
        let a_indices = c_indices.clone();
        let b_indices = c_indices
            .into_iter()
            .map(|idx| idx + cur_codeword.len() / 2)
            .collect::<Vec<_>>();

        let revealed_codewords_a = (0..self.num_colinearity_tests)
            .map(|s| cur_codeword[a_indices[s]].clone())
            .collect();
        let revealed_codewords_b = (0..self.num_colinearity_tests)
            .map(|s| cur_codeword[b_indices[s]].clone())
            .collect();
        let revealed_codewords_c = (0..self.num_colinearity_tests)
            .map(|s| next_codeword[c_indices[s]].clone())
            .collect();

        let serialized_cur_codeword = cur_codeword
            .into_iter()
            .map(|c| bincode::serialize(&c).expect("Serialization failed"))
            .collect::<Vec<_>>();
        let serialized_next_codeword = next_codeword
            .into_iter()
            .map(|c| bincode::serialize(&c).expect("Serialization failed"))
            .collect::<Vec<_>>();

        let revealed_paths_a = (0..self.num_colinearity_tests)
            .map(|s| Merkle::open(a_indices[s], &serialized_cur_codeword))
            .collect();
        let revealed_paths_b = (0..self.num_colinearity_tests)
            .map(|s| Merkle::open(b_indices[s], &serialized_cur_codeword))
            .collect();
        let revealed_paths_c = (0..self.num_colinearity_tests)
            .map(|s| Merkle::open(c_indices[s], &serialized_next_codeword))
            .collect();

        FriQueryLayer {
            a: (revealed_codewords_a, revealed_paths_a),
            b: (revealed_codewords_b, revealed_paths_b),
            c: (revealed_codewords_c, revealed_paths_c),
        }
    }

    pub fn verify(&self, proof: &FriProof<F>, polynomial_values: &mut Vec<(usize, F)>) -> bool {
        let mut proof_stream = FiatShamirTransformer::new();

        let mut omega = self.omega.clone();
        let mut offset = self.offset.clone();

        // extract all roots and alphas
        let mut alphas: Vec<F> = Vec::new();
        for r in &proof.merkle_roots {
            proof_stream.push(&vec![r.to_vec()]);
            alphas.push(F::sample(&proof_stream.prover_fiat_shamir(32)));
        }

        // extract last codeword
        let serialized_last_codeword = proof
            .last_codeword
            .clone()
            .into_iter()
            .map(|c| bincode::serialize(&c).expect("Serialization failed"))
            .collect::<Vec<_>>();
        proof_stream.push(&serialized_last_codeword);

        // check if it matches the given root
        if **proof.merkle_roots.last().unwrap() != Merkle::commit(&serialized_last_codeword) {
            return false;
        }

        // check if it is low degree
        let degree = (proof.last_codeword.len() / self.expansion_factor) - 1;
        let mut last_omega = omega.clone();
        let mut last_offset = offset.clone();
        for r in 0..(self.num_rounds() - 1) {
            last_omega = last_omega.clone() * last_omega;
            last_offset = last_offset.clone() * last_offset;
        }

        // assert that last_omega has the right order
        assert!(
            last_omega.inverse() == last_omega.pow(proof.last_codeword.len() - 1),
            "omega does not have right order"
        );

        // compute interpolant
        let last_domain = (0..proof.last_codeword.len())
            .into_iter()
            .map(|i| last_offset.clone() * (last_omega.pow(i)))
            .collect::<Vec<_>>();
        let poly = Polynomial::<F>::interpolate(&last_domain, &proof.last_codeword);

        for i in 0..proof.last_codeword.len() {
            assert!(
                poly.eval(&last_domain[i]) == proof.last_codeword[i],
                "re-evaluated codeword does not match original!"
            );
        }

        if poly.degree() > degree.try_into().unwrap() {
            return false;
        }

        let top_level_indices = sample_indices(
            &proof_stream.prover_fiat_shamir(32),
            self.domain_length >> 1,
            self.domain_length >> (self.num_rounds() - 1),
            self.num_colinearity_tests,
        );

        for r in 0..(self.num_rounds() - 1) {
            let c_indices = top_level_indices
                .iter()
                .map(|i| i % (self.domain_length >> (r + 1)))
                .collect::<Vec<_>>();
            let a_indices = c_indices.clone();
            let b_indices = a_indices
                .iter()
                .map(|i| i + (self.domain_length >> (r + 1)))
                .collect::<Vec<_>>();

            let mut aa = Vec::<F>::new();
            let mut bb = Vec::<F>::new();
            let mut cc = Vec::<F>::new();
            for s in 0..self.num_colinearity_tests {
                let ay: F = proof.revealed_layers[r].a.0[s].clone();
                let by: F = proof.revealed_layers[r].b.0[s].clone();
                let cy: F = proof.revealed_layers[r].c.0[s].clone();
                aa.push(ay.clone());
                bb.push(by.clone());
                cc.push(cy.clone());

                if r == 0 {
                    polynomial_values.push((a_indices[s], ay.clone()));
                    polynomial_values.push((b_indices[s], by.clone()));
                    //polynomial_values.push((c_indices[s], cy.clone()));
                }

                // colinearity check
                let ax = &offset.mul_ref(&(omega.pow(a_indices[s])));
                let bx = &offset.mul_ref(&(omega.pow(b_indices[s])));
                let cx = &alphas[r];
                let p = Polynomial::<F>::interpolate(
                    &[ax.clone(), bx.clone(), cx.clone()],
                    &[ay.clone(), by.clone(), cy.clone()],
                );
                if p.degree() > 1 {
                    return false;
                }
            }

            // verify authentication paths
            for i in 0..self.num_colinearity_tests {
                let path_aa = proof.revealed_layers[r].a.1[i].clone();
                if !Merkle::verify(
                    &proof.merkle_roots[r],
                    a_indices[i],
                    &path_aa,
                    &bincode::serialize(&aa[i]).expect("Serialization failed"),
                ) {
                    println!("merkle authentication path verification fails for aa");
                    return false;
                }
                let path_bb = proof.revealed_layers[r].b.1[i].clone();
                if !Merkle::verify(
                    &proof.merkle_roots[r],
                    b_indices[i],
                    &path_bb,
                    &bincode::serialize(&bb[i]).expect("Serialization failed"),
                ) {
                    println!("merkle authentication path verification fails for bb");
                    return false;
                }
                let path_cc = proof.revealed_layers[r].c.1[i].clone();
                if !Merkle::verify(
                    &proof.merkle_roots[r + 1],
                    c_indices[i],
                    &path_cc,
                    &bincode::serialize(&cc[i]).expect("Serialization failed"),
                ) {
                    println!("merkle authentication path verification fails for cc");
                    return false;
                }
            }

            omega = omega.pow(2);
            offset = offset.pow(2);
        }

        true
    }
}

define_myzkp_modulus_type!(M128, "270497897142230380135924736767050121217");
define_myzkp_modulus_type!(M64, "18446744069414584321");
define_extension_field!(
    Ip3,
    FiniteFieldElement<M64>,
    Polynomial {
        coef: vec![
            FiniteFieldElement::<M64>::one(),
            FiniteFieldElement::<M64>::zero() - FiniteFieldElement::<M64>::one(),
            FiniteFieldElement::<M64>::zero(),
            FiniteFieldElement::<M64>::one(),
        ],
    }
);

pub fn get_nth_root_of_m128(n: &BigInt) -> FiniteFieldElement<M128> {
    let one = BigInt::one();
    let shift_119 = one.clone().shl(119u32);
    // let expected_p = &one.clone() + &BigInt::from(407u32) * &shift_119;

    //assert_eq!(self.p, expected_p, "Field is not the expected one");

    // Check that n <= 2^119 and n is a power of two
    assert!(
        n <= &shift_119 && (n & (n - &one)).is_zero(),
        "Field does not have nth root of unity where n > 2^119 or not power of two."
    );

    let mut root = FiniteFieldElement::<M128>::from_value(
        BigInt::from_str("85408008396924667383611388730472331217").unwrap(),
    );
    let mut order = shift_119.clone();

    while order != *n {
        root = &root * &root; // or root = &root * &root;
        order = order.shr(1u32);
    }

    root
}

pub fn get_nth_root_of_m64(n: &BigInt) -> ExtendedFieldElement<M64, Ip3> {
    let one = BigInt::one();
    let shift_119 = one.clone().shl(32u32);
    // let expected_p = &one.clone() + &BigInt::from(407u32) * &shift_119;

    //assert_eq!(self.p, expected_p, "Field is not the expected one");

    // Check that n <= 2^119 and n is a power of two
    assert!(
        n <= &shift_119 && (n & (n - &one)).is_zero(),
        "Field does not have nth root of unity where n > 2^119 or not power of two."
    );

    let mut root = ExtendedFieldElement::<M64, Ip3>::from_value(
        BigInt::from_str("1753635133440165772").unwrap(),
    );
    let mut order = shift_119.clone();

    while order != *n {
        root = root.mul_ref(&root); // or root = &root * &root;
        order = order.shr(1u32);
    }

    root
}

pub fn compute_log_codeword_length(initial_codeword_length: usize) -> usize {
    let mut codeword_length = initial_codeword_length;
    let mut log_codeword_length = 0;

    while codeword_length > 1 {
        codeword_length /= 2;
        log_codeword_length += 1;
    }

    log_codeword_length
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::modules::algebra::field::FiniteFieldElement;

    // 1 + 407 * (1 << 119)

    #[test]
    fn test_fri_field() {
        let degree: usize = 63;
        let expansion_factor: usize = 4;
        let num_colinearity_tests: usize = 17;
        let initial_codeword_length = (degree + 1) * expansion_factor;

        let log_codeword_length = compute_log_codeword_length(initial_codeword_length);
        assert_eq!(1 << log_codeword_length, initial_codeword_length);

        let generator = FiniteFieldElement::<M128>::from_value(
            BigInt::from_str("85408008396924667383611388730472331217").unwrap(),
        );
        let omega = get_nth_root_of_m128(&BigInt::from(initial_codeword_length));
        assert!(omega.pow(1 << log_codeword_length).is_one());
        assert!(!omega.pow(1 << (log_codeword_length - 1)).is_one());

        let fri = FRI {
            offset: generator,
            omega: omega.clone(),
            domain_length: initial_codeword_length,
            expansion_factor: expansion_factor,
            num_colinearity_tests: num_colinearity_tests,
        };

        let polynomial = Polynomial {
            coef: (0..degree + 1)
                .map(|i| FiniteFieldElement::<M128>::from_value(i))
                .collect(),
        };

        let domain: Vec<_> = (0..initial_codeword_length).map(|i| omega.pow(i)).collect();
        let mut codeword = polynomial.eval_domain(&domain);
        let proof = fri.prove(&codeword);

        let mut points = Vec::new();
        let result = fri.verify(&proof, &mut points);
        assert!(result);
        for (x, y) in &points {
            assert_eq!(polynomial.eval(&omega.pow(*x)), y.clone());
        }

        for i in 0..(degree / 3) {
            codeword[i] = FiniteFieldElement::<M128>::one();
        }
        let proof = fri.prove(&codeword);
        let mut points_second = Vec::new();
        let result_second = fri.verify(&proof, &mut points_second);
        assert!(!result_second);
    }

    #[test]
    fn test_fri_efield() {
        let degree: usize = 63;
        let expansion_factor: usize = 16;
        let num_colinearity_tests: usize = 17;
        let initial_codeword_length = (degree + 1) * expansion_factor;

        let log_codeword_length = compute_log_codeword_length(initial_codeword_length);
        assert_eq!(1 << log_codeword_length, initial_codeword_length);

        let generator =
            ExtendedFieldElement::<M64, Ip3>::from_value(BigInt::from_str("7").unwrap());
        let omega = get_nth_root_of_m64(&BigInt::from(initial_codeword_length));
        assert!(omega.pow(1 << log_codeword_length).is_one());
        //assert!(!omega.pow(1 << (log_codeword_length - 1)).is_one());

        let fri = FRI {
            offset: generator,
            omega: omega.clone(),
            domain_length: initial_codeword_length,
            expansion_factor: expansion_factor,
            num_colinearity_tests: num_colinearity_tests,
        };

        let polynomial = Polynomial {
            coef: (0..degree + 1)
                .map(|i| ExtendedFieldElement::<M64, Ip3>::from_value(i))
                .collect(),
        };

        let domain: Vec<_> = (0..initial_codeword_length).map(|i| omega.pow(i)).collect();
        let mut codeword = polynomial.eval_domain(&domain);
        let proof = fri.prove(&codeword);

        let mut points = Vec::new();
        let result = fri.verify(&proof, &mut points);
        assert!(result);
        for (x, y) in &points {
            assert_eq!(polynomial.eval(&omega.pow(*x)), y.clone());
        }

        for i in 0..(degree / 3) {
            codeword[i] = ExtendedFieldElement::<M64, Ip3>::one();
        }
        let proof = fri.prove(&codeword);
        let mut points_second = Vec::new();
        let result_second = fri.verify(&proof, &mut points_second);
        assert!(!result_second);
    }
}
