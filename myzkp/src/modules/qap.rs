use num_traits::One;
use num_traits::Zero;

use crate::modules::field::Field;
use crate::modules::polynomial::Polynomial;
use crate::modules::r1cs::{dot, R1CS};

#[derive(Debug, Clone)]
pub struct QAP<'a, F: Field> {
    pub r1cs: &'a R1CS<F>,
    pub t: Polynomial<F>,
}

impl<'a, F: Field> QAP<'a, F> {
    fn new(r1cs: &'a R1CS<F>) -> Self {
        QAP {
            r1cs: r1cs,
            t: Polynomial::<F>::from_monomials(
                &(1..=r1cs.d).map(|i| F::from_value(i)).collect::<Vec<F>>(),
            ),
        }
    }

    fn generate_polynomials(&self, a: &Vec<F>) -> (Polynomial<F>, Polynomial<F>, Polynomial<F>) {
        let left_dot_products = self
            .r1cs
            .left
            .iter()
            .map(|v| dot(&v, &a))
            .collect::<Vec<F>>();
        let right_dot_products = self
            .r1cs
            .right
            .iter()
            .map(|v| dot(&v, &a))
            .collect::<Vec<F>>();
        let out_dot_products = self
            .r1cs
            .out
            .iter()
            .map(|v| dot(&v, &a))
            .collect::<Vec<F>>();

        let x = (1..=self.r1cs.m)
            .map(|i| F::from_value(i))
            .collect::<Vec<F>>();
        let left_interpolated_polynomial = Polynomial::<F>::interpolate(&x, &left_dot_products);
        let right_interpolated_polynomial = Polynomial::<F>::interpolate(&x, &right_dot_products);
        let out_interpolated_polynomial = Polynomial::<F>::interpolate(&x, &out_dot_products);
        (
            left_interpolated_polynomial,
            right_interpolated_polynomial,
            out_interpolated_polynomial,
        )
    }
}
