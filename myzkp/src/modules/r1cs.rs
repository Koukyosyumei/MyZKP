use num_traits::One;
use num_traits::Zero;

use crate::modules::field::Field;

pub fn dot<F: Field>(a: &Vec<F>, b: &Vec<F>) -> F {
    let mut result = F::zero();
    for (a_i, b_i) in a.iter().zip(b.iter()) {
        result = result + a_i.clone() * b_i.clone();
    }
    result
}

#[derive(Debug, Clone)]
pub struct R1CS<F: Field> {
    pub left: Vec<Vec<F>>,
    pub right: Vec<Vec<F>>,
    pub out: Vec<Vec<F>>,
    pub m: usize,
    pub d: usize,
}

impl<F: Field> R1CS<F> {
    pub fn new(left: Vec<Vec<F>>, right: Vec<Vec<F>>, out: Vec<Vec<F>>) -> Self {
        let d = left.len();
        let m = if d == 0 { 0 } else { left[0].len() };
        R1CS {
            left,
            right,
            out,
            m,
            d,
        }
    }

    pub fn is_satisfied(&self, a: &Vec<F>) -> bool {
        let zero = F::zero();
        self.left
            .iter()
            .zip(self.right.iter())
            .zip(self.out.iter())
            .all(|((l, r), o)| dot(&l, &a) * dot(&r, &a) - dot(&o, &a) == zero)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::modules::field::{FiniteFieldElement, ModEIP197};
    use crate::modules::ring::Ring;

    type F = FiniteFieldElement<ModEIP197>;

    #[test]
    fn test_r1cs_single_multiplication() {
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
        let a_prime = vec![
            F::one(),
            F::from_value(3690),
            F::from_value(82),
            F::from_value(46),
        ];
        let r1cs = R1CS::new(left, right, out);
        assert!(r1cs.is_satisfied(&a));
        assert!(!r1cs.is_satisfied(&a_prime));
    }

    #[test]
    fn test_r1cs_addition_with_constant() {
        // z = x * y + 3
        // (1, z, x, y) = (1, 9, 2, 3)
        let left = vec![vec![F::zero(), F::zero(), F::one(), F::zero()]];
        let right = vec![vec![F::zero(), F::zero(), F::zero(), F::one()]];
        let out = vec![vec![F::from_value(-3_i32), F::one(), F::zero(), F::zero()]];
        let a = vec![
            F::one(),
            F::from_value(9),
            F::from_value(2),
            F::from_value(3),
        ];
        let a_prime = vec![
            F::one(),
            F::from_value(9),
            F::from_value(2),
            F::from_value(4),
        ];
        let r1cs = R1CS::new(left, right, out);
        assert!(r1cs.is_satisfied(&a));
        assert!(!r1cs.is_satisfied(&a_prime));
    }

    #[test]
    fn test_r1cs_multiplication_with_constant() {
        // z = 3x^2 + y
        // (1, z, x, y) = (1, 17, 2, 5)
        let left = vec![vec![F::zero(), F::zero(), F::from_value(3_u32), F::zero()]];
        let right = vec![vec![F::zero(), F::zero(), F::one(), F::zero()]];
        let out = vec![vec![F::zero(), F::one(), F::zero(), F::from_value(-1_i32)]];
        let a = vec![
            F::one(),
            F::from_value(17),
            F::from_value(2),
            F::from_value(5),
        ];
        let a_prime = vec![
            F::one(),
            F::from_value(17),
            F::from_value(2),
            F::from_value(7),
        ];
        let r1cs = R1CS::new(left, right, out);
        assert!(r1cs.is_satisfied(&a));
        assert!(!r1cs.is_satisfied(&a_prime));
    }

    #[test]
    fn test_r1cs_multiple_constraints() {
        // v = a * b
        // u = c * d
        // r = v * u
        // (1, r, a, b, c, d, v, u) = (1, 210, 2, 3, 5, 7, 6, 35)
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

        let a = vec![
            F::one(),
            F::from_value(210),
            F::from_value(2),
            F::from_value(3),
            F::from_value(5),
            F::from_value(7),
            F::from_value(6),
            F::from_value(35),
        ];
        let a_prime = vec![
            F::one(),
            F::from_value(211),
            F::from_value(2),
            F::from_value(3),
            F::from_value(5),
            F::from_value(7),
            F::from_value(6),
            F::from_value(35),
        ];
        let r1cs = R1CS::new(left, right, out);
        assert!(r1cs.is_satisfied(&a));
        assert!(!r1cs.is_satisfied(&a_prime));
    }
}
