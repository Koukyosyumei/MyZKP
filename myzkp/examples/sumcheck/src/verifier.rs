
use cudarc;

use myzkp::modules::algebra::fiat_shamir::FiatShamirTransformer;
use myzkp::modules::algebra::field::Field;
use myzkp::modules::algebra::mpolynomials::MPolynomial;
use myzkp::modules::algebra::polynomial::Polynomial;
use myzkp::modules::algebra::ring::Ring;

use crate::utils::{F};

pub struct SumCheckVerifier;

impl SumCheckVerifier {
    pub fn verify(
        max_degree: usize,
        polynomial_factors: &[MPolynomial<F>],
        claimed_sum: F,
        proof: &Vec<u8>,
    ) -> Result<bool, Box<dyn std::error::Error>> {
        let mut transcript = FiatShamirTransformer::deserialize(&proof);

        let num_variables = polynomial_factors.iter().map(|p| p.get_num_vars()).max().unwrap_or(0);
        let polynomial_product = polynomial_factors.iter().skip(1).fold(polynomial_factors[0].clone(), |acc, p| &acc * p);
        let s_eval_point: Vec<F> = (0..(max_degree + 1)).map(|d| F::from_value(d)).collect();

        let recovered_max_degree = transcript.pull();
        assert_eq!(
            vec![bincode::serialize(&max_degree).expect("Serialization failed")],
            recovered_max_degree
        );
        let recovered_num_factors = bincode::deserialize::<usize>(&transcript.pull()[0])?;
        assert_eq!(recovered_num_factors, polynomial_factors.len());
        let recovered_num_variables = bincode::deserialize::<usize>(&transcript.pull()[0])?;
        assert_eq!(recovered_num_variables, num_variables);
        for i in 0..recovered_num_factors {
            assert_eq!(
                vec![bincode::serialize(&polynomial_factors[i]).expect("Serialization failed")],
                transcript.pull()
            );
        }

        let mut s_polys = vec![];
        let mut challenges: Vec<F> = vec![];
        for i in 0..num_variables {
            let mut s_evals = vec![];
            for d in 0..(max_degree + 1) {
                let tmp_bytes = transcript.pull();
                let se = bincode::deserialize::<F>(&tmp_bytes[0])?;
                s_evals.push(se);
            }

            let s_poly = Polynomial::interpolate(&s_eval_point, &s_evals);
            s_polys.push(s_poly);

            let challenge = F::sample(&transcript.verifier_fiat_shamir(32));
            challenges.push(challenge);

            if i == 0 {
                assert_eq!(s_evals[0].add_ref(&s_evals[1]), claimed_sum);
            } else {
                assert_eq!(
                    s_evals[0].add_ref(&s_evals[1]),
                    s_polys[i - 1].eval(&challenges[i - 1]).sanitize()
                );
            }
        }
        assert_eq!(
            polynomial_product.evaluate(&challenges),
            s_polys[num_variables - 1]
                .eval(&challenges[num_variables - 1])
                .sanitize()
        );

        Ok(true)
    }
}
