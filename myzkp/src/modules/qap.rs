use num_traits::One;
use num_traits::Zero;

use crate::modules::field::Field;
use crate::modules::polynomial::Polynomial;
use crate::modules::r1cs::R1CS;

pub struct QAP<F: Field> {
    pub m: usize,
    pub d: usize,
    pub ell_i_vec: Vec<Polynomial<F>>,
    pub r_i_vec: Vec<Polynomial<F>>,
    pub o_i_vec: Vec<Polynomial<F>>,
    pub t: Polynomial<F>,
}

impl<F: Field> QAP<F> {
    pub fn from_r1cs(r1cs: &R1CS<F>) -> QAP<F> {
        let x: Vec<F> = (1..=r1cs.m).map(F::from_value).collect();

        let interpolate = |rows: &[Vec<F>]| {
            (0..r1cs.d)
                .map(|i| {
                    let column: Vec<F> = rows.iter().map(|row| row[i].clone()).collect();
                    Polynomial::<F>::interpolate(&x, &column)
                })
                .collect::<Vec<Polynomial<F>>>()
        };

        let ell_i_vec = interpolate(&r1cs.left);
        let r_i_vec = interpolate(&r1cs.right);
        let o_i_vec = interpolate(&r1cs.out);

        let t = Polynomial::<F>::from_monomials(&x);

        QAP {
            m: r1cs.m,
            d: r1cs.d,
            ell_i_vec,
            r_i_vec,
            o_i_vec,
            t,
        }
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
        let r1cs = R1CS::new(left, right, out);
        let qap = QAP::from_r1cs(&r1cs);

        let ground_truth_ell_i_vec = vec![
            Polynomial::<F>::zero(),
            Polynomial::<F>::zero(),
            Polynomial::<F>::one(),
            Polynomial::<F>::zero(),
        ];
        let ground_truth_r_i_vec = vec![
            Polynomial::<F>::zero(),
            Polynomial::<F>::zero(),
            Polynomial::<F>::zero(),
            Polynomial::<F>::one(),
        ];
        let ground_truth_o_i_vec = vec![
            Polynomial::<F>::zero(),
            Polynomial::<F>::one(),
            Polynomial::<F>::zero(),
            Polynomial::<F>::zero(),
        ];
        for i in 0..(qap.d) {
            assert_eq!(qap.ell_i_vec[i], ground_truth_ell_i_vec[i]);
            assert_eq!(qap.r_i_vec[i], ground_truth_r_i_vec[i]);
            assert_eq!(qap.o_i_vec[i], ground_truth_o_i_vec[i]);
        }
    }

    #[test]
    fn test_qap_multi_multiplication() {
        let left = vec![
            vec![
                F::zero(),
                F::zero(),
                F::one(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
            ],
            vec![
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::one(),
                F::zero(),
                F::zero(),
                F::zero(),
            ],
            vec![
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::one(),
                F::zero(),
            ],
        ];

        let right = vec![
            vec![
                F::zero(),
                F::zero(),
                F::zero(),
                F::one(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
            ],
            vec![
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::one(),
                F::zero(),
                F::zero(),
            ],
            vec![
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::one(),
            ],
        ];

        let out = vec![
            vec![
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::one(),
                F::zero(),
            ],
            vec![
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::one(),
            ],
            vec![
                F::zero(),
                F::one(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
            ],
        ];
        let r1cs = R1CS::new(left, right, out);
        let qap = QAP::from_r1cs(&r1cs);

        /*
        for i in 0..(qap.d) {
            println!("ell - {}: {}", i, qap.ell_i_vec[i]);
        }

        assert!(false);
        */
    }
}
