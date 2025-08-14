use std::fmt;

use crate::modules::algebra::field::Field;
use crate::modules::algebra::polynomial::Polynomial;

#[derive(Debug)]
pub enum SplitFoldError {
    PointsLenMismatch { expected: usize, found: usize },
    CoefsNotPowerOfTwo { found: usize },
}

impl fmt::Display for SplitFoldError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SplitFoldError::PointsLenMismatch { expected, found } => {
                write!(f, "points.len() must be {}, but got {}", expected, found)
            }
            SplitFoldError::CoefsNotPowerOfTwo { found } => {
                write!(f, "coefs.len() must be a power of two, but got {}", found)
            }
        }
    }
}

fn int_log2(n: usize) -> u32 {
    debug_assert!(n > 0);
    usize::BITS - 1 - n.leading_zeros()
}

pub fn tensor_product<F: Field>(a: &Vec<F>, b: &Vec<F>) -> Vec<F> {
    // a * b^T
    let ab = a
        .iter()
        .map(|v| b.iter().map(|u| v.mul_ref(u)).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    (0..b.len())
        .map(|i| ab.iter().map(|vs| vs[i].clone()).collect::<Vec<_>>())
        .flatten()
        .collect::<Vec<_>>()
}

pub fn split_and_fold<F: Field>(
    coef: &Vec<F>,
    rhos: &Vec<F>,
) -> Result<Vec<Polynomial<F>>, SplitFoldError> {
    let n = coef.len();
    if n.count_ones() != 1 {
        return Err(SplitFoldError::CoefsNotPowerOfTwo { found: n });
    }

    let log2_n = int_log2(n);
    if rhos.len() != log2_n as usize {
        return Err(SplitFoldError::PointsLenMismatch {
            expected: log2_n as usize,
            found: rhos.len(),
        });
    }

    let mut f = Polynomial::<F> { coef: coef.clone() };
    let mut fs = vec![f.clone()];
    for i in 1..(log2_n + 1) {
        let f_e = Polynomial::<F> {
            coef: f
                .coef
                .iter()
                .enumerate()
                .map(|(i, x)| if i % 2 == 0 { x.clone() } else { F::zero() })
                .collect(),
        };
        let f_o = Polynomial::<F> {
            coef: f
                .coef
                .iter()
                .skip(1)
                .enumerate()
                .map(|(i, x)| if i % 2 == 0 { x.clone() } else { F::zero() })
                .collect(),
        };
        let f_i = f_e + f_o * rhos[i as usize - 1].clone();

        f = Polynomial::<F> {
            coef: f_i
                .coef
                .iter()
                .enumerate()
                .filter(|(i, _)| i % 2 == 0)
                .map(|(_, x)| x.clone())
                .collect(),
        };
        fs.push(f.clone());
    }

    Ok(fs)
}

pub fn vanila_verify<F: Field>(
    rhos: &Vec<F>,
    mu: &F,
    polys: &Vec<Polynomial<F>>,
    beta: &F,
) -> bool {
    let log2_n = rhos.len();
    let es = polys
        .iter()
        .take(log2_n)
        .map(|f| f.eval(beta))
        .collect::<Vec<_>>();
    let es_neg = polys
        .iter()
        .take(log2_n)
        .map(|f| f.eval(&(F::zero() - beta)))
        .collect::<Vec<_>>();
    let mut es_hat = polys
        .iter()
        .skip(1)
        .take(log2_n - 1)
        .map(|f| f.eval(&beta.pow(2)))
        .collect::<Vec<_>>();
    es_hat.push(mu.clone());

    let two = F::from_value(2);
    (0..log2_n).all(|j| {
        two.mul_ref(beta).mul_ref(&es_hat[j])
            == beta.mul_ref(&es[j].add_ref(&es_neg[j]))
                + rhos[j].mul_ref(&es[j].sub_ref(&es_neg[j]))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    use num_traits::{One, Zero};

    use crate::modules::algebra::field::{FiniteFieldElement, ModEIP197};
    use crate::modules::algebra::ring::Ring;

    type F = FiniteFieldElement<ModEIP197>;

    #[test]
    fn test_tensor_product_1() {
        let a = vec![F::from_value(2), F::from_value(3)];
        let b = vec![F::from_value(4), F::from_value(5), F::from_value(6)];
        let c = tensor_product(&a, &b);
        let c_groundtruth = vec![
            F::from_value(8),
            F::from_value(12),
            F::from_value(10),
            F::from_value(15),
            F::from_value(12),
            F::from_value(18),
        ];
        assert_eq!(c.len(), c_groundtruth.len());
        for (u, v) in c.iter().zip(c_groundtruth.iter()) {
            assert_eq!(u, v);
        }
    }

    #[test]
    fn test_tensor_product_2() {
        let x1 = vec![F::from_value(1), F::from_value(2)];
        let x2 = vec![F::from_value(1), F::from_value(3)];
        let x3 = vec![F::from_value(1), F::from_value(4)];
        let c = tensor_product(&tensor_product(&x1, &x2), &x3);
        let c_groundtruth = vec![
            F::from_value(1),
            F::from_value(2),
            F::from_value(3),
            F::from_value(6),
            F::from_value(4),
            F::from_value(8),
            F::from_value(12),
            F::from_value(24),
        ];
        assert_eq!(c.len(), c_groundtruth.len());
        for (u, v) in c.iter().zip(c_groundtruth.iter()) {
            assert_eq!(u, v);
        }
    }

    #[test]
    fn test_vanila_verifier() {
        let coef = (0..8).map(|i| F::from_value(i + 1)).collect::<Vec<_>>();
        let rhos = vec![F::from_value(2), F::from_value(3), F::from_value(4)];
        let c = tensor_product(
            &tensor_product(
                &vec![F::one(), rhos[0].clone()],
                &vec![F::one(), rhos[1].clone()],
            ),
            &vec![F::one(), rhos[2].clone()],
        );

        let mut mu = F::zero();
        for (u, v) in coef.iter().zip(&c) {
            mu += u.mul_ref(v);
        }

        let fs = split_and_fold(&coef, &rhos).unwrap();
        assert!(vanila_verify(&rhos, &mu, &fs, &F::from_value(1234)));
        assert!(!vanila_verify(
            &rhos,
            &(mu + F::one()),
            &fs,
            &F::from_value(1234)
        ));
    }
}
