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

#[cfg(test)]
mod tests {
    use super::*;

    use crate::modules::field::{FiniteFieldElement, ModEIP197};
    use crate::modules::ring::Ring;

    type F = FiniteFieldElement<ModEIP197>;

    #[test]
    fn test_qap_single_multiplication() {
        // z = x * y
        // (1, z, x, y) = (1, 3690, 82, 45)
        let left = vec![vec![F::zero(), F::zero(), F::one(), F::zero()]];
        let right = vec![vec![F::zero(), F::zero(), F::zero(), F::one()]];
        let out = vec![vec![F::zero(), F::one(), F::zero(), F::zero()]];
        let a = vec![
            F::one(),
            F::from_value(3690),
            F::from_value(82),
            F::from_value(45),
        ];
        let r1cs = R1CS::new(left, right, out);
        let qap = QAP::new(&r1cs);
        let (_, _, _) = qap.generate_polynomials(&a);
    }
}
