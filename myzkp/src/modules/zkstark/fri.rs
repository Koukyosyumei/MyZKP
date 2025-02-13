use crate::modules::algebra::field::Field;
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::merkle::merkle::Merkle;
use crate::modules::zkstark::fiat_shamir::FiatShamirTransformer;

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

    fn sample_indices(&self) -> Vec<usize> {
        todo!()
    }

    pub fn prove(&self, codeword: &Codeword<F>, proof_stream: &mut FiatShamirTransformer) {
        let codewords = self.commit(codeword, proof_stream);
        let top_level_indices = self.sample_indices();
        let mut indices = top_level_indices.clone();

        for i in 0..(codewords.len() - 1) {
            indices = indices
                .into_iter()
                .map(|idx| idx % (codewords[i].len() / 2))
                .collect::<Vec<_>>()
                .clone();
            self.query(&codewords[i], &codewords[i + 1], &indices, proof_stream);
        }
    }

    pub fn query(
        &self,
        cur_codeword: &Codeword<F>,
        next_codeword: &Codeword<F>,
        c_indices: &Vec<usize>,
        proof_stream: &FiatShamirTransformer,
    ) {
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
            let root = Merkle::commit(
                &codeword
                    .clone()
                    .into_iter()
                    .map(|c| c.serialize())
                    .collect::<Vec<_>>(),
            );
            proof_stream.push(&root);

            if r == self.num_rounds() - 1 {
                break;
            }

            let alpha = F::one();
            codewords.push(codeword.clone());

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

        // proof_stream.push(codeword);
        codewords.push(codeword);
        codewords
    }
}
