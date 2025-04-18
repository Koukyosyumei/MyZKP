use std::collections::HashMap;
use std::fmt::Debug;
use std::ops::{Add, Mul, Neg, Sub};

use num_traits::Zero;

use crate::modules::algebra::field::Field;
use crate::modules::algebra::polynomial::Polynomial;

#[derive(Clone, Debug)]
pub struct MPolynomial<F: Field> {
    pub dictionary: HashMap<Vec<usize>, F>,
}

impl<F: Field> MPolynomial<F> {
    pub fn new(dictionary: HashMap<Vec<usize>, F>) -> Self {
        // Filter out zero coefficients
        let filtered: HashMap<Vec<usize>, F> = dictionary
            .into_iter()
            .filter(|(_, v)| !v.is_zero())
            .collect();

        MPolynomial {
            dictionary: filtered,
        }
    }

    pub fn zero() -> Self {
        MPolynomial {
            dictionary: HashMap::new(),
        }
    }

    pub fn is_zero(&self) -> bool {
        self.dictionary.is_empty() || self.dictionary.values().all(|v| v.is_zero())
    }

    pub fn constant(element: F) -> Self {
        if element.is_zero() {
            return MPolynomial::zero();
        }

        let mut dictionary = HashMap::new();
        dictionary.insert(vec![0], element);
        MPolynomial { dictionary }
    }

    pub fn variables(num_variables: usize) -> Vec<Self> {
        let mut variables = Vec::with_capacity(num_variables);

        for i in 0..num_variables {
            let mut exponent = vec![0; num_variables];
            exponent[i] = 1;

            let mut dictionary = HashMap::new();
            dictionary.insert(exponent, F::one());

            variables.push(MPolynomial { dictionary });
        }

        variables
    }

    // Implementation of the Python __xor__ method (exponentiation)
    pub fn pow(&self, exponent: usize) -> Self {
        if self.is_zero() {
            return MPolynomial::zero();
        }

        // Get the number of variables
        let num_variables = self.dictionary.keys().next().map_or(1, |k| k.len());

        // Create identity polynomial (1)
        let mut acc_dict = HashMap::new();
        acc_dict.insert(vec![0; num_variables], F::one());
        let mut acc = MPolynomial {
            dictionary: acc_dict,
        };

        // Binary exponentiation algorithm
        for bit in format!("{:b}", exponent).chars() {
            acc = &acc * &acc;

            if bit == '1' {
                acc = &acc * self;
            }
        }

        acc
    }

    pub fn evaluate(&self, point: &[F]) -> F {
        let mut acc = F::zero();

        for (k, v) in &self.dictionary {
            let mut prod = v.clone();

            for i in 0..k.len() {
                // Skip if i is out of bounds in point
                if i >= point.len() {
                    continue;
                }

                // Calculate point[i]^k[i]
                let mut term = F::one();
                for _ in 0..k[i] {
                    term = term * point[i].clone();
                }

                prod = prod * term;
            }

            acc = acc + prod;
        }

        acc
    }

    pub fn evaluate_symbolic(&self, point: &[Polynomial<F>]) -> Polynomial<F> {
        let mut acc = Polynomial::zero();

        for (k, v) in &self.dictionary {
            let mut prod = Polynomial {
                coef: [v.clone()].to_vec(),
            };

            for i in 0..k.len() {
                prod = prod * (point[i].pow(k[i]));
            }

            acc = acc + prod;
        }

        acc
    }

    pub fn lift(polynomial: &Polynomial<F>, variable_index: usize) -> Self {
        if polynomial.is_zero() {
            return MPolynomial::zero();
        }

        let coefficients = polynomial.coef.clone();
        // Create variables x_0, x_1, ..., x_{variable_index}
        let variables = MPolynomial::variables(variable_index + 1);

        // Get the last variable (x_{variable_index})
        let x = &variables[variable_index];

        let mut acc = MPolynomial::zero();
        for (i, coef) in coefficients.iter().enumerate() {
            let constant_term = MPolynomial::constant(coef.clone());
            // Multiply by coef * x^i
            let term = &constant_term * &x.pow(i);
            acc = &acc + &term;
        }

        acc
    }
}

// Addition implementation
impl<F: Field> Add for &MPolynomial<F> {
    type Output = MPolynomial<F>;

    fn add(self, other: &MPolynomial<F>) -> MPolynomial<F> {
        let mut dictionary = HashMap::new();

        // Find the maximum number of variables
        let num_variables = self
            .dictionary
            .keys()
            .map(|k| k.len())
            .chain(other.dictionary.keys().map(|k| k.len()))
            .max()
            .unwrap_or(0);

        // Process self's terms
        for (k, v) in &self.dictionary {
            let mut pad = k.clone();
            pad.resize(num_variables, 0);
            dictionary.insert(pad, v.clone());
        }

        // Process other's terms
        for (k, v) in &other.dictionary {
            let mut pad = k.clone();
            pad.resize(num_variables, 0);

            if let Some(existing) = dictionary.get_mut(&pad) {
                *existing = existing.clone() + v.clone();
                if existing.is_zero() {
                    dictionary.remove(&pad);
                }
            } else {
                dictionary.insert(pad, v.clone());
            }
        }

        MPolynomial { dictionary }
    }
}

impl<F: Field> Add for MPolynomial<F> {
    type Output = MPolynomial<F>;

    fn add(self, other: MPolynomial<F>) -> MPolynomial<F> {
        &self + &other
    }
}

// Multiplication implementation
impl<F: Field> Mul for &MPolynomial<F> {
    type Output = MPolynomial<F>;

    fn mul(self, other: &MPolynomial<F>) -> MPolynomial<F> {
        let mut dictionary: HashMap<Vec<usize>, F> = HashMap::new();

        // Find the maximum number of variables
        let num_variables = self
            .dictionary
            .keys()
            .map(|k| k.len())
            .chain(other.dictionary.keys().map(|k| k.len()))
            .max()
            .unwrap_or(0);

        for (k0, v0) in &self.dictionary {
            for (k1, v1) in &other.dictionary {
                let mut exponent = vec![0; num_variables];

                // Add exponents from k0
                for (k, &v) in k0.iter().enumerate() {
                    if k < exponent.len() {
                        exponent[k] += v;
                    }
                }

                // Add exponents from k1
                for (k, &v) in k1.iter().enumerate() {
                    if k < exponent.len() {
                        exponent[k] += v;
                    }
                }

                let product = v0.clone() * v1.clone();

                if let Some(existing) = dictionary.get_mut(&exponent) {
                    *existing = existing.clone() + product;
                    if existing.is_zero() {
                        dictionary.remove(&exponent);
                    }
                } else {
                    dictionary.insert(exponent, product);
                }
            }
        }

        MPolynomial { dictionary }
    }
}

impl<F: Field> Mul for MPolynomial<F> {
    type Output = MPolynomial<F>;

    fn mul(self, other: MPolynomial<F>) -> MPolynomial<F> {
        &self * &other
    }
}

// Negation implementation
impl<F: Field> Neg for &MPolynomial<F> {
    type Output = MPolynomial<F>;

    fn neg(self) -> MPolynomial<F> {
        let dictionary = self
            .dictionary
            .iter()
            .map(|(k, v)| (k.clone(), -v.clone()))
            .collect();

        MPolynomial { dictionary }
    }
}

impl<F: Field> Neg for MPolynomial<F> {
    type Output = MPolynomial<F>;

    fn neg(self) -> MPolynomial<F> {
        -&self
    }
}

// Subtraction implementation
impl<F: Field> Sub for &MPolynomial<F> {
    type Output = MPolynomial<F>;

    fn sub(self, other: &MPolynomial<F>) -> MPolynomial<F> {
        self + &(-other)
    }
}

impl<F: Field> Sub for MPolynomial<F> {
    type Output = MPolynomial<F>;

    fn sub(self, other: MPolynomial<F>) -> MPolynomial<F> {
        &self - &other
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::modules::algebra::field::{FiniteFieldElement, ModEIP197};
    use crate::modules::algebra::ring::Ring;

    #[test]
    fn test_mpolynomial_addition() {
        // Create polynomials: 2x + 3y and 4x + 5z
        let mut dict1 = HashMap::new();
        dict1.insert(
            vec![1, 0, 0],
            FiniteFieldElement::<ModEIP197>::from_value(2),
        ); // 2x
        dict1.insert(
            vec![0, 1, 0],
            FiniteFieldElement::<ModEIP197>::from_value(3),
        ); // 3y
        let poly1 = MPolynomial::new(dict1);

        let mut dict2 = HashMap::new();
        dict2.insert(
            vec![1, 0, 0],
            FiniteFieldElement::<ModEIP197>::from_value(4),
        ); // 4x
        dict2.insert(
            vec![0, 0, 1],
            FiniteFieldElement::<ModEIP197>::from_value(5),
        ); // 5z
        let poly2 = MPolynomial::new(dict2);

        // Expected: 6x + 3y + 5z
        let result = &poly1 + &poly2;

        assert_eq!(result.dictionary.len(), 3);
        assert_eq!(
            result.dictionary[&vec![1, 0, 0]],
            FiniteFieldElement::<ModEIP197>::from_value(6)
        );
        assert_eq!(
            result.dictionary[&vec![0, 1, 0]],
            FiniteFieldElement::<ModEIP197>::from_value(3)
        );
        assert_eq!(
            result.dictionary[&vec![0, 0, 1]],
            FiniteFieldElement::<ModEIP197>::from_value(5)
        );
    }

    #[test]
    fn test_mpolynomial_addition_with_cancellation() {
        // Create polynomials: 2x + 3y and -2x + 5z
        let mut dict1 = HashMap::new();
        dict1.insert(
            vec![1, 0, 0],
            FiniteFieldElement::<ModEIP197>::from_value(2),
        );
        dict1.insert(
            vec![0, 1, 0],
            FiniteFieldElement::<ModEIP197>::from_value(3),
        );
        let poly1 = MPolynomial::new(dict1);

        let mut dict2 = HashMap::new();
        dict2.insert(
            vec![1, 0, 0],
            FiniteFieldElement::<ModEIP197>::from_value(-2),
        );
        dict2.insert(
            vec![0, 0, 1],
            FiniteFieldElement::<ModEIP197>::from_value(5),
        );
        let poly2 = MPolynomial::new(dict2);

        // Expected: 3y + 5z (x term cancels)
        let result = &poly1 + &poly2;

        assert_eq!(result.dictionary.len(), 2);
        assert!(!result.dictionary.contains_key(&vec![1, 0, 0]));
        assert_eq!(
            result.dictionary[&vec![0, 1, 0]],
            FiniteFieldElement::<ModEIP197>::from_value(3)
        );
        assert_eq!(
            result.dictionary[&vec![0, 0, 1]],
            FiniteFieldElement::<ModEIP197>::from_value(5)
        );
    }

    #[test]
    fn test_mpolynomial_multiplication() {
        // Create polynomials: 2x and 3y
        let mut dict1 = HashMap::new();
        dict1.insert(vec![1, 0], FiniteFieldElement::<ModEIP197>::from_value(2));
        let poly1 = MPolynomial::new(dict1);

        let mut dict2 = HashMap::new();
        dict2.insert(vec![0, 1], FiniteFieldElement::<ModEIP197>::from_value(3));
        let poly2 = MPolynomial::new(dict2);

        // Expected: 6xy
        let result = &poly1 * &poly2;

        assert_eq!(result.dictionary.len(), 1);
        assert_eq!(
            result.dictionary[&vec![1, 1]],
            FiniteFieldElement::<ModEIP197>::from_value(6)
        );
    }

    #[test]
    fn test_mpolynomial_multiplication_complex() {
        // Create polynomials: (x + 2) and (y + 3)
        let mut dict1 = HashMap::new();
        dict1.insert(vec![1, 0], FiniteFieldElement::<ModEIP197>::from_value(1)); // x
        dict1.insert(vec![0, 0], FiniteFieldElement::<ModEIP197>::from_value(2)); // 2
        let poly1 = MPolynomial::new(dict1);

        let mut dict2 = HashMap::new();
        dict2.insert(vec![0, 1], FiniteFieldElement::<ModEIP197>::from_value(1)); // y
        dict2.insert(vec![0, 0], FiniteFieldElement::<ModEIP197>::from_value(3)); // 3
        let poly2 = MPolynomial::new(dict2);

        // Expected: xy + 3x + 2y + 6
        let result = &poly1 * &poly2;

        assert_eq!(result.dictionary.len(), 4);
        assert_eq!(
            result.dictionary[&vec![1, 1]],
            FiniteFieldElement::<ModEIP197>::from_value(1)
        ); // xy
        assert_eq!(
            result.dictionary[&vec![1, 0]],
            FiniteFieldElement::<ModEIP197>::from_value(3)
        ); // 3x
        assert_eq!(
            result.dictionary[&vec![0, 1]],
            FiniteFieldElement::<ModEIP197>::from_value(2)
        ); // 2y
        assert_eq!(
            result.dictionary[&vec![0, 0]],
            FiniteFieldElement::<ModEIP197>::from_value(6)
        ); // 6
    }

    #[test]
    fn test_mpolynomial_evaluate() {
        // Create polynomial: 2x + 3y + 1
        let mut dict = HashMap::new();
        dict.insert(vec![1, 0], FiniteFieldElement::<ModEIP197>::from_value(2)); // 2x
        dict.insert(vec![0, 1], FiniteFieldElement::<ModEIP197>::from_value(3)); // 3y
        dict.insert(vec![0, 0], FiniteFieldElement::<ModEIP197>::from_value(1)); // 1
        let poly = MPolynomial::new(dict);

        // Evaluate at point (5, 7)
        let point = vec![
            FiniteFieldElement::<ModEIP197>::from_value(5),
            FiniteFieldElement::<ModEIP197>::from_value(7),
        ];
        let result = poly.evaluate(&point);

        // Expected: 2*5 + 3*7 + 1 = 10 + 21 + 1 = 32
        assert_eq!(result, FiniteFieldElement::<ModEIP197>::from_value(32));
    }

    #[test]
    fn test_evaluate_with_higher_degrees() {
        // Create polynomial: x^2 + 2xy + y^3
        let mut dict = HashMap::new();
        dict.insert(vec![2, 0], FiniteFieldElement::<ModEIP197>::from_value(1)); // x^2
        dict.insert(vec![1, 1], FiniteFieldElement::<ModEIP197>::from_value(2)); // 2xy
        dict.insert(vec![0, 3], FiniteFieldElement::<ModEIP197>::from_value(1)); // y^3
        let poly = MPolynomial::new(dict);

        // Evaluate at point (3, 2)
        let point = vec![
            FiniteFieldElement::<ModEIP197>::from_value(3),
            FiniteFieldElement::<ModEIP197>::from_value(2),
        ];
        let result = poly.evaluate(&point);

        // Expected: 3^2 + 2*3*2 + 2^3 = 9 + 12 + 8 = 29
        assert_eq!(result, FiniteFieldElement::<ModEIP197>::from_value(29));
    }

    #[test]
    fn test_evaluate_symbolic() {
        // Create polynomial: 2x + 3y
        let mut dict = HashMap::new();
        dict.insert(vec![1, 0], FiniteFieldElement::<ModEIP197>::from_value(2)); // 2x
        dict.insert(vec![0, 1], FiniteFieldElement::<ModEIP197>::from_value(3)); // 3y
        let poly = MPolynomial::new(dict);

        // Symbolic points: x = (t + 1), y = (t^2)
        let point = vec![
            Polynomial {
                coef: vec![
                    FiniteFieldElement::<ModEIP197>::from_value(1),
                    FiniteFieldElement::<ModEIP197>::from_value(1),
                ],
            }, // t + 1
            Polynomial {
                coef: vec![
                    FiniteFieldElement::<ModEIP197>::from_value(0),
                    FiniteFieldElement::<ModEIP197>::from_value(0),
                    FiniteFieldElement::<ModEIP197>::from_value(1),
                ],
            }, // t^2
        ];

        let result = poly.evaluate_symbolic(&point);

        // Expected: 2(t + 1) + 3t^2 = 2 + 2t + 3t^2
        assert_eq!(result.coef.len(), 3);
        assert_eq!(
            result.coef[0],
            FiniteFieldElement::<ModEIP197>::from_value(2)
        );
        assert_eq!(
            result.coef[1],
            FiniteFieldElement::<ModEIP197>::from_value(2)
        );
        assert_eq!(
            result.coef[2],
            FiniteFieldElement::<ModEIP197>::from_value(3)
        );
    }

    #[test]
    fn test_evaluate_symbolic_constant() {
        // Create constant polynomial: 5
        let mut dict = HashMap::new();
        dict.insert(vec![0, 0], FiniteFieldElement::<ModEIP197>::from_value(5));
        let poly = MPolynomial::new(dict);

        // Symbolic points: x = (t^2), y = (2t + 3)
        let point = vec![
            Polynomial {
                coef: vec![
                    FiniteFieldElement::<ModEIP197>::from_value(0),
                    FiniteFieldElement::<ModEIP197>::from_value(0),
                    FiniteFieldElement::<ModEIP197>::from_value(1),
                ],
            }, // t^2
            Polynomial {
                coef: vec![
                    FiniteFieldElement::<ModEIP197>::from_value(3),
                    FiniteFieldElement::<ModEIP197>::from_value(2),
                ],
            }, // 2t + 3
        ];

        let result = poly.evaluate_symbolic(&point);

        // Expected: 5 (constant)
        assert_eq!(result.coef.len(), 1);
        assert_eq!(
            result.coef[0],
            FiniteFieldElement::<ModEIP197>::from_value(5)
        );
    }

    #[test]
    fn test_lift() {
        // Create a univariate polynomial: 1 + 2t + 3t^2
        let uni_poly = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(1),
                FiniteFieldElement::<ModEIP197>::from_value(2),
                FiniteFieldElement::<ModEIP197>::from_value(3),
            ],
        };

        // Lift to a multivariate polynomial with variable index 2
        // This should create a polynomial in terms of x_2: 1 + 2x_2 + 3x_2^2
        let multi_poly = MPolynomial::lift(&uni_poly, 2);

        assert_eq!(multi_poly.dictionary.len(), 3);
        assert_eq!(
            multi_poly.dictionary[&vec![0, 0, 0]],
            FiniteFieldElement::<ModEIP197>::from_value(1)
        ); // Constant term
        assert_eq!(
            multi_poly.dictionary[&vec![0, 0, 1]],
            FiniteFieldElement::<ModEIP197>::from_value(2)
        ); // x_2 term
        assert_eq!(
            multi_poly.dictionary[&vec![0, 0, 2]],
            FiniteFieldElement::<ModEIP197>::from_value(3)
        ); // x_2^2 term
    }
}
