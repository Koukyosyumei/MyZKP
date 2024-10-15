use num_bigint::ToBigInt;
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

// Assuming FiniteFieldElement is already implemented with necessary traits like Add, Sub, Mul, Div.
use crate::modules::field::Field;
use crate::modules::field::FiniteFieldElement;

/// Polynomial struct representing a polynomial over FiniteFieldElement.
#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial {
    pub poly: Vec<FiniteFieldElement>,
    pub var: String, // Variable for representation, default to 'x'
}

impl Polynomial {
    /// Creates a polynomial representing the variable `x`.
    pub fn x() -> Self {
        Polynomial {
            poly: vec![
                FiniteFieldElement::zero(None),
                FiniteFieldElement::one(None),
            ],
            var: "x".to_string(),
        }
    }

    /// Removes trailing zeroes from a polynomial's coefficients.
    fn trim_trailing_zeros(poly: Vec<FiniteFieldElement>) -> Vec<FiniteFieldElement> {
        let mut trimmed = poly;
        while trimmed.last() == Some(&FiniteFieldElement::zero(None)) {
            trimmed.pop();
        }
        trimmed
    }

    /// Returns the degree of the polynomial.
    pub fn degree(&self) -> isize {
        let trimmed = Self::trim_trailing_zeros(self.poly.clone());
        if trimmed.is_empty() {
            -1
        } else {
            (trimmed.len() - 1) as isize
        }
    }

    /// Returns the nth coefficient.
    pub fn nth_coefficient(&self, n: usize) -> FiniteFieldElement {
        if n > self.degree() as usize {
            FiniteFieldElement::zero(None)
        } else {
            self.poly[n].clone()
        }
    }

    /// Scalar multiplication of a polynomial.
    pub fn scalar_mul(&self, scalar: &FiniteFieldElement) -> Polynomial {
        Polynomial {
            poly: self
                .poly
                .iter()
                .map(|x| x.clone() * scalar.clone())
                .collect(),
            var: self.var.clone(),
        }
    }

    /// Evaluate the polynomial at a given point.
    pub fn eval(&self, point: &FiniteFieldElement) -> FiniteFieldElement {
        let mut result = FiniteFieldElement::zero(None);
        for coef in self.poly.iter().rev() {
            result = result * point.clone() + coef.clone();
        }
        result
    }

    /// Lagrange interpolation to compute polynomials.
    pub fn interpolate(
        x_values: &[FiniteFieldElement],
        y_values: &[FiniteFieldElement],
    ) -> Polynomial {
        let mut lagrange_polys = vec![];
        let numerators = Polynomial::from_monomials(x_values);

        for j in 0..x_values.len() {
            // \Pi_{j \neq i} (x_j - x_i)
            let mut denominator = FiniteFieldElement::one(None);
            for i in 0..x_values.len() {
                if i != j {
                    denominator = denominator * (x_values[j].clone() - x_values[i].clone());
                }
            }
            let cur_poly = numerators
                .clone()
                .div(Polynomial::from_monomials(&[x_values[j].clone()]).scalar_mul(&denominator));
            lagrange_polys.push(cur_poly);
        }

        let mut result = Polynomial {
            poly: vec![],
            var: "x".to_string(),
        };
        for (j, lagrange_poly) in lagrange_polys.iter().enumerate() {
            result = result + lagrange_poly.scalar_mul(&y_values[j]);
        }
        result
    }

    /// Helper to create polynomial from a single monomial.
    pub fn from_monomials(x_values: &[FiniteFieldElement]) -> Polynomial {
        let mut poly = Polynomial {
            poly: vec![FiniteFieldElement::one(None)],
            var: "x".to_string(),
        };
        for x in x_values {
            poly = poly.mul(Polynomial {
                poly: vec![
                    FiniteFieldElement::zero(None) - x.clone(),
                    FiniteFieldElement::one(None),
                ],
                var: "x".to_string(),
            });
        }
        poly
    }
}

impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut terms = Vec::new();

        for (i, coeff) in self.poly.iter().enumerate() {
            if coeff != &FiniteFieldElement::zero(None) {
                let term = if i == 0 {
                    // Constant term
                    format!("{}", coeff)
                } else if i == 1 {
                    // Linear term (e.g., 3x)
                    if coeff == &FiniteFieldElement::one(None) {
                        format!("{}", self.var)
                    } else {
                        format!("{}{}", coeff, self.var)
                    }
                } else {
                    // Higher degree terms (e.g., 3x^2)
                    if coeff == &FiniteFieldElement::one(None) {
                        format!("{}^{}", self.var, i)
                    } else {
                        format!("{}{}^{}", coeff, self.var, i)
                    }
                };
                terms.push(term);
            }
        }

        // If there are no non-zero terms, return "0"
        if terms.is_empty() {
            write!(f, "0")
        } else {
            // Join the terms with " + " and print the result
            write!(f, "{}", terms.join(" + "))
        }
    }
}

// Arithmetic operations implementation for Polynomial.
impl Add for Polynomial {
    type Output = Self;

    fn add(self, other: Self) -> Polynomial {
        let max_len = std::cmp::max(self.poly.len(), other.poly.len());
        let mut result = Vec::with_capacity(max_len);

        let zero = FiniteFieldElement::zero(None);

        for i in 0..max_len {
            let a = self.poly.get(i).unwrap_or(&zero);
            let b = other.poly.get(i).unwrap_or(&zero);
            result.push(a.clone() + b.clone());
        }
        Polynomial {
            poly: Self::trim_trailing_zeros(result),
            var: self.var.clone(),
        }
    }
}

impl Sub for Polynomial {
    type Output = Self;

    fn sub(self, other: Self) -> Polynomial {
        let max_len = std::cmp::max(self.poly.len(), other.poly.len());
        let mut result = Vec::with_capacity(max_len);

        let zero = FiniteFieldElement::zero(None);

        for i in 0..max_len {
            let a = self.poly.get(i).unwrap_or(&zero);
            let b = other.poly.get(i).unwrap_or(&zero);
            result.push(a.clone() - b.clone());
        }
        Polynomial {
            poly: Self::trim_trailing_zeros(result),
            var: self.var.clone(),
        }
    }
}

impl Neg for Polynomial {
    type Output = Self;

    /// Negation of a polynomial.
    fn neg(self) -> Polynomial {
        Polynomial {
            poly: self.poly.iter().map(|x| -x.clone()).collect(),
            var: self.var.clone(),
        }
    }
}

impl Mul for Polynomial {
    type Output = Self;

    /// Multiplication of two polynomials.
    fn mul(self, other: Polynomial) -> Polynomial {
        let mut result = vec![
            FiniteFieldElement::zero(None);
            self.degree() as usize + other.degree() as usize + 1
        ];

        for (i, a) in self.poly.iter().enumerate() {
            for (j, b) in other.poly.iter().enumerate() {
                result[i + j] = result[i + j].clone() + (a.clone() * b.clone());
            }
        }
        Polynomial {
            poly: Self::trim_trailing_zeros(result),
            var: self.var.clone(),
        }
    }
}

impl Div for Polynomial {
    type Output = Self;

    /// Division of two polynomials, returns quotient.
    fn div(self, other: Polynomial) -> Polynomial {
        let mut remainder_coeffs = Self::trim_trailing_zeros(self.poly.clone());
        let divisor_coeffs = Self::trim_trailing_zeros(other.poly.clone());
        let divisor_lead_inv = divisor_coeffs.last().unwrap().inverse();

        let mut quotient = vec![
            FiniteFieldElement::zero(None);
            self.degree() as usize - other.degree() as usize + 1
        ];

        let mut i = 0_i32;
        while remainder_coeffs.len() >= divisor_coeffs.len() && i < 20 {
            let lead_term = remainder_coeffs.last().unwrap().clone() * divisor_lead_inv.clone();
            let deg_diff = remainder_coeffs.len() - divisor_coeffs.len();
            quotient[deg_diff] = lead_term.clone();

            for i in 0..divisor_coeffs.len() {
                remainder_coeffs[deg_diff + i] = remainder_coeffs[deg_diff + i].clone()
                    - (lead_term.clone() * divisor_coeffs[i].clone());
            }
            remainder_coeffs = Self::trim_trailing_zeros(remainder_coeffs);
            i += 1;
        }

        Polynomial {
            poly: Self::trim_trailing_zeros(quotient),
            var: self.var.clone(),
        }

        /*
        (
            Polynomial {
                poly: Self::trim_trailing_zeros(quotient),
                var: self.var.clone(),
            },
            Polynomial {
                poly: remainder,
                var: self.var.clone(),
            },
        )
        */
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polynomial_addition() {
        let poly1 = Polynomial {
            poly: vec![
                FiniteFieldElement::from(1_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(2_i32.to_bigint().unwrap()),
            ],
            var: "x".to_string(),
        };
        let poly2 = Polynomial {
            poly: vec![
                FiniteFieldElement::from(2_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(3_i32.to_bigint().unwrap()),
            ],
            var: "x".to_string(),
        };
        let expected = Polynomial {
            poly: vec![
                FiniteFieldElement::from(3_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(5_i32.to_bigint().unwrap()),
            ],
            var: "x".to_string(),
        };
        assert_eq!(poly1 + poly2, expected);
    }

    #[test]
    fn test_polynomial_subtraction() {
        let poly1 = Polynomial {
            poly: vec![
                FiniteFieldElement::from(4_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(5_i32.to_bigint().unwrap()),
            ],
            var: "x".to_string(),
        };
        let poly2 = Polynomial {
            poly: vec![
                FiniteFieldElement::from(1_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(3_i32.to_bigint().unwrap()),
            ],
            var: "x".to_string(),
        };
        let expected = Polynomial {
            poly: vec![
                FiniteFieldElement::from(3_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(2_i32.to_bigint().unwrap()),
            ],
            var: "x".to_string(),
        };
        assert_eq!(poly1 - poly2, expected);
    }

    #[test]
    fn test_polynomial_negation() {
        let poly = Polynomial {
            poly: vec![
                FiniteFieldElement::from(3_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(4_i32.to_bigint().unwrap()),
            ],
            var: "x".to_string(),
        };
        let expected = Polynomial {
            poly: vec![
                FiniteFieldElement::from(-3_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(-4_i32.to_bigint().unwrap()),
            ],
            var: "x".to_string(),
        };
        assert_eq!(-poly, expected);
    }

    #[test]
    fn test_polynomial_multiplication() {
        let poly1 = Polynomial {
            poly: vec![
                FiniteFieldElement::from(1_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(2_i32.to_bigint().unwrap()),
            ], // 1 + 2x
            var: "x".to_string(),
        };
        let poly2 = Polynomial {
            poly: vec![
                FiniteFieldElement::from(2_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(3_i32.to_bigint().unwrap()),
            ], // 2 + 3x
            var: "x".to_string(),
        };
        let expected = Polynomial {
            poly: vec![
                FiniteFieldElement::from(2_i32.to_bigint().unwrap()), // constant term
                FiniteFieldElement::from(7_i32.to_bigint().unwrap()), // x term
                FiniteFieldElement::from(6_i32.to_bigint().unwrap()), // x^2 term
            ],
            var: "x".to_string(),
        };
        assert_eq!(poly1 * poly2, expected);
    }

    #[test]
    fn test_polynomial_scalar_multiplication() {
        let poly = Polynomial {
            poly: vec![
                FiniteFieldElement::from(1_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(2_i32.to_bigint().unwrap()),
            ], // 1 + 2x
            var: "x".to_string(),
        };
        let scalar = FiniteFieldElement::from(3_i32.to_bigint().unwrap());
        let expected = Polynomial {
            poly: vec![
                FiniteFieldElement::from(3_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(6_i32.to_bigint().unwrap()),
            ], // 3 + 6x
            var: "x".to_string(),
        };
        assert_eq!(poly.scalar_mul(&scalar), expected);
    }

    #[test]
    fn test_polynomial_division() {
        let poly1 = Polynomial {
            poly: vec![
                FiniteFieldElement::from(2_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(3_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(1_i32.to_bigint().unwrap()),
            ], // 2 + 3x + x^2
            var: "x".to_string(),
        };
        let poly2 = Polynomial {
            poly: vec![
                FiniteFieldElement::from(1_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(1_i32.to_bigint().unwrap()),
            ], // 1 + x
            var: "x".to_string(),
        };
        let quotient = poly1 / poly2;
        let expected_quotient = Polynomial {
            poly: vec![
                FiniteFieldElement::from(2_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(1_i32.to_bigint().unwrap()),
            ], // 2 + x
            var: "x".to_string(),
        };
        //let expected_remainder = Polynomial {
        //    poly: vec![FiniteFieldElement::from(1_i32.to_bigint().unwrap())], // remainder is 1
        //    var: "x".to_string(),
        //};
        assert_eq!(quotient, expected_quotient);
        // assert_eq!(remainder, expected_remainder);
    }

    #[test]
    fn test_polynomial_evaluation() {
        let poly = Polynomial {
            poly: vec![
                FiniteFieldElement::from(2_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(3_i32.to_bigint().unwrap()),
            ], // 2 + 3x
            var: "x".to_string(),
        };
        let point = FiniteFieldElement::from(2_i32.to_bigint().unwrap());
        let expected = FiniteFieldElement::from(8_i32.to_bigint().unwrap()); // 2 + 3 * 2 = 8
        assert_eq!(poly.eval(&point), expected);
    }

    #[test]
    fn test_polynomial_degree() {
        let poly = Polynomial {
            poly: vec![
                FiniteFieldElement::from(1_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(0_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(3_i32.to_bigint().unwrap()),
            ], // 1 + 0x + 3x^2
            var: "x".to_string(),
        };
        assert_eq!(poly.degree(), 2); // The degree should be 2
    }

    #[test]
    fn test_polynomial_lagrange_interpolation() {
        let x_values = vec![
            FiniteFieldElement::from(1_i32.to_bigint().unwrap()),
            FiniteFieldElement::from(2_i32.to_bigint().unwrap()),
            FiniteFieldElement::from(3_i32.to_bigint().unwrap()),
        ];
        let y_values = vec![
            FiniteFieldElement::from(0_i32.to_bigint().unwrap()),
            FiniteFieldElement::from(3_i32.to_bigint().unwrap()),
            FiniteFieldElement::from(8_i32.to_bigint().unwrap()),
        ];
        let result = Polynomial::interpolate(&x_values, &y_values);
        let expected = Polynomial {
            poly: vec![
                FiniteFieldElement::from(-1_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(0_i32.to_bigint().unwrap()),
                FiniteFieldElement::from(1_i32.to_bigint().unwrap()),
            ], // x^2 - 1
            var: "x".to_string(),
        };
        assert_eq!(result, expected);
    }

    #[test]
    fn test_polynomial_from_monomials() {
        let points = vec![
            FiniteFieldElement::from(2_i32.to_bigint().unwrap()),
            FiniteFieldElement::from(3_i32.to_bigint().unwrap()),
        ];
        let result = Polynomial::from_monomials(&points);
        // (x - 2) * (x - 3) = x^2 - 5x + 6
        let expected = Polynomial {
            poly: vec![
                FiniteFieldElement::from(6_i32.to_bigint().unwrap()), // constant term
                FiniteFieldElement::from(-5_i32.to_bigint().unwrap()), // x term
                FiniteFieldElement::from(1_i32.to_bigint().unwrap()), // x^2 term
            ],
            var: "x".to_string(),
        };
        assert_eq!(result, expected);
    }
}
