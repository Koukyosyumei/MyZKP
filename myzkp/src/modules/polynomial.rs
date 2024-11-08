use num_traits::{One, Zero};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Rem, Sub};

use crate::modules::curve::{EllipticCurve, EllipticCurvePoint};
use crate::modules::field::Field;

/// Polynomial struct representing a polynomial over Field.
#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial<F: Field> {
    pub coef: Vec<F>,
}

impl<F: Field> Polynomial<F> {
    /// Creates a polynomial representing the variable `x`.
    pub fn x() -> Self {
        Polynomial {
            coef: vec![F::zero(), F::one()],
        }
    }

    /// Removes trailing zeroes from a polynomial's coefficients.
    fn trim_trailing_zeros(coef: Vec<F>) -> Vec<F> {
        let mut trimmed = coef;
        while trimmed.last() == Some(&F::zero()) {
            trimmed.pop();
        }
        trimmed
    }

    pub fn reduce(&self) -> Self {
        Polynomial {
            coef: Self::trim_trailing_zeros(self.coef.clone()),
        }
    }

    /// Returns the degree of the polynomial.
    pub fn degree(&self) -> isize {
        let trimmed = Self::trim_trailing_zeros(self.coef.clone());
        if trimmed.is_empty() {
            -1
        } else {
            (trimmed.len() - 1) as isize
        }
    }

    /// Returns the nth coefficient.
    pub fn nth_coefficient(&self, n: usize) -> F {
        if n > self.degree() as usize {
            F::zero()
        } else {
            self.coef[n].clone()
        }
    }

    /// Evaluate the polynomial at a given point.
    pub fn eval(&self, point: &F) -> F {
        let mut result = F::zero();
        for coef in self.coef.iter().rev() {
            result = result.mul_ref(&point) + coef;
        }
        result
    }

    pub fn eval_with_powers(&self, powers: &[F]) -> F {
        let mut result = F::one();
        for (i, coef) in self.coef.iter().enumerate() {
            result = result * powers[i].pow(coef.clone().get_value());
        }
        result
    }

    pub fn eval_with_powers_on_curve<E: EllipticCurve>(
        &self,
        powers: &[EllipticCurvePoint<F, E>],
    ) -> EllipticCurvePoint<F, E> {
        let mut result = EllipticCurvePoint::point_at_infinity();
        for (i, coef) in self.coef.iter().enumerate() {
            result = result + powers[i].clone() * coef.clone().get_value();
        }
        result
    }

    /// Lagrange interpolation to compute polynomials.
    pub fn interpolate(x_values: &[F], y_values: &[F]) -> Polynomial<F> {
        let mut lagrange_polys = vec![];
        let numerators = Polynomial::from_monomials(x_values);

        for j in 0..x_values.len() {
            let mut denominator = F::one();
            for i in 0..x_values.len() {
                if i != j {
                    denominator = denominator * (x_values[j].clone() - x_values[i].clone());
                }
            }
            let cur_poly = numerators
                .clone()
                .div(Polynomial::from_monomials(&[x_values[j].clone()]) * denominator);
            lagrange_polys.push(cur_poly);
        }

        let mut result = Polynomial { coef: vec![] };
        for (j, lagrange_poly) in lagrange_polys.iter().enumerate() {
            result = result + lagrange_poly.clone() * y_values[j].clone();
        }
        result
    }

    /// Helper to create polynomial from a single monomial.
    pub fn from_monomials(x_values: &[F]) -> Polynomial<F> {
        let mut poly = Polynomial {
            coef: vec![F::one()],
        };
        for x in x_values {
            poly = poly.mul(Polynomial {
                coef: vec![F::zero() - x.clone(), F::one()],
            });
        }
        poly
    }

    fn add_ref<'a, 'b>(&self, other: &'b Polynomial<F>) -> Polynomial<F> {
        let max_len = std::cmp::max(self.coef.len(), other.coef.len());
        let mut result = Vec::with_capacity(max_len);

        let zero = F::zero();

        for i in 0..max_len {
            let a = self.coef.get(i).unwrap_or(&zero);
            let b = other.coef.get(i).unwrap_or(&zero);
            result.push(a.clone() + b.clone());
        }
        Polynomial {
            coef: Self::trim_trailing_zeros(result),
        }
    }

    fn mul_ref<'a, 'b>(&self, other: &'b Polynomial<F>) -> Polynomial<F> {
        if self.is_zero() || other.is_zero() {
            return Polynomial::<F>::zero();
        }
        let mut result = vec![F::zero(); (self.degree() + other.degree() + 1) as usize];

        for (i, a) in self.coef.iter().enumerate() {
            for (j, b) in other.coef.iter().enumerate() {
                result[i + j] = result[i + j].clone() + (a.clone() * b.clone());
            }
        }
        Polynomial {
            coef: Polynomial::<F>::trim_trailing_zeros(result),
        }
    }
}

impl<F: Field> fmt::Display for Polynomial<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut terms = Vec::new();

        for (i, coeff) in self.coef.iter().enumerate() {
            if coeff != &F::zero() {
                let term = if i == 0 {
                    // Constant term
                    format!("{}", coeff)
                } else if i == 1 {
                    // Linear term (e.g., 3x)
                    if coeff == &F::one() {
                        format!("{}", "x")
                    } else {
                        format!("{}{}", coeff, "x")
                    }
                } else {
                    // Higher degree terms (e.g., 3x^2)
                    if coeff == &F::one() {
                        format!("{}^{}", "x", i)
                    } else {
                        format!("{}{}^{}", coeff, "x", i)
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

impl<F: Field> Zero for Polynomial<F> {
    fn zero() -> Self {
        Polynomial {
            coef: vec![F::zero()],
        }
    }

    fn is_zero(&self) -> bool {
        self.degree() == -1
    }
}

impl<F: Field> One for Polynomial<F> {
    fn one() -> Self {
        Polynomial {
            coef: vec![F::one()],
        }
    }
}

// Arithmetic operations implementation for Polynomial.
impl<F: Field> Neg for Polynomial<F> {
    type Output = Self;

    fn neg(self) -> Polynomial<F> {
        Polynomial {
            coef: self.coef.iter().map(|x| -x.clone()).collect(),
        }
    }
}

impl<F: Field> Add for Polynomial<F> {
    type Output = Self;

    fn add(self, other: Self) -> Polynomial<F> {
        self.add_ref(&other)
    }
}

impl<'a, 'b, F: Field> Add<&'b Polynomial<F>> for &'a Polynomial<F> {
    type Output = Polynomial<F>;

    fn add(self, other: &'b Polynomial<F>) -> Polynomial<F> {
        self.add_ref(other)
    }
}

impl<F: Field> Sub for Polynomial<F> {
    type Output = Self;

    fn sub(self, other: Self) -> Polynomial<F> {
        self.add_ref(&-other)
    }
}

impl<'a, 'b, F: Field> Sub<&'b Polynomial<F>> for &'a Polynomial<F> {
    type Output = Polynomial<F>;

    fn sub(self, other: &'b Polynomial<F>) -> Polynomial<F> {
        self.add_ref(&-other.clone())
    }
}

impl<F: Field> Mul<Polynomial<F>> for Polynomial<F> {
    type Output = Self;

    fn mul(self, other: Polynomial<F>) -> Polynomial<F> {
        self.mul_ref(&other)
    }
}

impl<'a, 'b, F: Field> Mul<&'b Polynomial<F>> for &'a Polynomial<F> {
    type Output = Polynomial<F>;

    fn mul(self, other: &'b Polynomial<F>) -> Polynomial<F> {
        self.mul_ref(other)
    }
}

impl<F: Field> Mul<F> for Polynomial<F> {
    type Output = Self;

    fn mul(self, scalar: F) -> Polynomial<F> {
        Polynomial {
            coef: self.coef.iter().map(|x| x.mul_ref(&scalar)).collect(),
        }
    }
}

impl<F: Field> Div for Polynomial<F> {
    type Output = Self;

    /// Division of two polynomials, returns quotient.
    fn div(self, other: Polynomial<F>) -> Polynomial<F> {
        if self.degree() < other.degree() {
            return self.clone();
        }

        let mut remainder_coeffs = Self::trim_trailing_zeros(self.coef.clone());
        let divisor_coeffs = Self::trim_trailing_zeros(other.coef.clone());
        let divisor_lead_inv = divisor_coeffs.last().unwrap().inverse();

        let mut quotient = vec![F::zero(); self.degree() as usize - other.degree() as usize + 1];

        while remainder_coeffs.len() >= divisor_coeffs.len() {
            let lead_term = remainder_coeffs.last().unwrap().clone() * divisor_lead_inv.clone();
            let deg_diff = remainder_coeffs.len() - divisor_coeffs.len();
            quotient[deg_diff] = lead_term.clone();

            for i in 0..divisor_coeffs.len() {
                remainder_coeffs[deg_diff + i] = remainder_coeffs[deg_diff + i].clone()
                    - (lead_term.mul_ref(&divisor_coeffs[i]));
            }
            remainder_coeffs = Self::trim_trailing_zeros(remainder_coeffs);
        }

        Polynomial {
            coef: Self::trim_trailing_zeros(quotient),
        }
    }
}

impl<F: Field> Rem for Polynomial<F> {
    type Output = Self;

    /// Division of two polynomials, returns quotient.
    fn rem(self, other: Polynomial<F>) -> Polynomial<F> {
        if self.degree() < other.degree() {
            return self.clone();
        }

        let mut remainder_coeffs = Self::trim_trailing_zeros(self.coef.clone());
        let divisor_coeffs = Self::trim_trailing_zeros(other.coef.clone());
        let divisor_lead_inv = divisor_coeffs.last().unwrap().inverse();

        let mut quotient = vec![F::zero(); self.degree() as usize - other.degree() as usize + 1];

        while remainder_coeffs.len() >= divisor_coeffs.len() {
            let lead_term = remainder_coeffs.last().unwrap().mul_ref(&divisor_lead_inv);
            let deg_diff = remainder_coeffs.len() - divisor_coeffs.len();
            quotient[deg_diff] = lead_term.clone();

            for i in 0..divisor_coeffs.len() {
                remainder_coeffs[deg_diff + i] = remainder_coeffs[deg_diff + i].clone()
                    - (lead_term.mul_ref(&divisor_coeffs[i]));
            }
            remainder_coeffs = Self::trim_trailing_zeros(remainder_coeffs);
        }

        Polynomial {
            coef: remainder_coeffs,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::modules::field::{FiniteFieldElement, ModEIP197};
    use crate::modules::ring::Ring;

    #[test]
    fn test_polynomial_addition() {
        let poly1 = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(1_i32),
                FiniteFieldElement::<ModEIP197>::from_value(2_i32),
            ],
        };
        let poly2 = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(2_i32),
                FiniteFieldElement::<ModEIP197>::from_value(3_i32),
            ],
        };
        let expected = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(3_i32),
                FiniteFieldElement::<ModEIP197>::from_value(5_i32),
            ],
        };
        assert_eq!(poly1 + poly2, expected);
    }

    #[test]
    fn test_polynomial_subtraction() {
        let poly1 = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(4_i32),
                FiniteFieldElement::<ModEIP197>::from_value(5_i32),
            ],
        };
        let poly2 = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(1_i32),
                FiniteFieldElement::<ModEIP197>::from_value(3_i32),
            ],
        };
        let expected = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(3_i32),
                FiniteFieldElement::<ModEIP197>::from_value(2_i32),
            ],
        };
        assert_eq!(poly1 - poly2, expected);
    }

    #[test]
    fn test_polynomial_negation() {
        let poly = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(3_i32),
                FiniteFieldElement::<ModEIP197>::from_value(4_i32),
            ],
        };
        let expected = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(-3_i32),
                FiniteFieldElement::<ModEIP197>::from_value(-4_i32),
            ],
        };
        assert_eq!(-poly, expected);
    }

    #[test]
    fn test_polynomial_multiplication() {
        let poly1 = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(1_i32),
                FiniteFieldElement::<ModEIP197>::from_value(2_i32),
            ], // 1 + 2x
        };
        let poly2 = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(2_i32),
                FiniteFieldElement::<ModEIP197>::from_value(3_i32),
            ], // 2 + 3x
        };
        let expected = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(2_i32), // constant term
                FiniteFieldElement::<ModEIP197>::from_value(7_i32), // x term
                FiniteFieldElement::<ModEIP197>::from_value(6_i32), // x^2 term
            ],
        };
        assert_eq!(poly1 * poly2, expected);
    }

    #[test]
    fn test_polynomial_scalar_multiplication() {
        let poly = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(1_i32),
                FiniteFieldElement::<ModEIP197>::from_value(2_i32),
            ], // 1 + 2x
        };
        let scalar = FiniteFieldElement::<ModEIP197>::from_value(3_i32);
        let expected = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(3_i32),
                FiniteFieldElement::<ModEIP197>::from_value(6_i32),
            ], // 3 + 6x
        };
        assert_eq!(poly * scalar, expected);
    }

    #[test]
    fn test_polynomial_division() {
        let poly1 = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(3_i32),
                FiniteFieldElement::<ModEIP197>::from_value(3_i32),
                FiniteFieldElement::<ModEIP197>::from_value(1_i32),
            ], // 3 + 3x + x^2
        };
        let poly2 = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(1_i32),
                FiniteFieldElement::<ModEIP197>::from_value(1_i32),
            ], // 1 + x
        };
        let quotient = poly1.clone() / poly2.clone();
        let expected_quotient = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(2_i32),
                FiniteFieldElement::<ModEIP197>::from_value(1_i32),
            ], // 2 + x
        };
        let remainder = poly1.clone() % poly2.clone();
        let expected_remainder = Polynomial {
            coef: vec![FiniteFieldElement::<ModEIP197>::from_value(1_i32)], // remainder is 1
        };
        assert_eq!(quotient, expected_quotient);
        assert_eq!(remainder, expected_remainder);
    }

    #[test]
    fn test_polynomial_evaluation() {
        let poly = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(2_i32),
                FiniteFieldElement::<ModEIP197>::from_value(3_i32),
            ], // 2 + 3x
        };
        let point = FiniteFieldElement::<ModEIP197>::from_value(2_i32);
        let expected = FiniteFieldElement::<ModEIP197>::from_value(8_i32); // 2 + 3 * 2 = 8
        assert_eq!(poly.eval(&point), expected);
    }

    #[test]
    fn test_polynomial_degree() {
        let poly = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(1_i32),
                FiniteFieldElement::<ModEIP197>::from_value(0_i32),
                FiniteFieldElement::<ModEIP197>::from_value(3_i32),
            ], // 1 + 0x + 3x^2
        };
        assert_eq!(poly.degree(), 2); // The degree should be 2
    }

    #[test]
    fn test_polynomial_lagrange_interpolation() {
        let x_values = vec![
            FiniteFieldElement::<ModEIP197>::from_value(1_i32),
            FiniteFieldElement::<ModEIP197>::from_value(2_i32),
            FiniteFieldElement::<ModEIP197>::from_value(3_i32),
        ];
        let y_values = vec![
            FiniteFieldElement::<ModEIP197>::from_value(0_i32),
            FiniteFieldElement::<ModEIP197>::from_value(3_i32),
            FiniteFieldElement::<ModEIP197>::from_value(8_i32),
        ];
        let result = Polynomial::interpolate(&x_values, &y_values);
        let expected = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(-1_i32).sanitize(),
                FiniteFieldElement::<ModEIP197>::from_value(0_i32),
                FiniteFieldElement::<ModEIP197>::from_value(1_i32),
            ], // x^2 - 1
        };
        assert_eq!(result, expected);
    }

    #[test]
    fn test_polynomial_from_monomials() {
        let points = vec![
            FiniteFieldElement::<ModEIP197>::from_value(2_i32),
            FiniteFieldElement::<ModEIP197>::from_value(3_i32),
        ];
        let result = Polynomial::from_monomials(&points);
        // (x - 2) * (x - 3) = x^2 - 5x + 6
        let expected = Polynomial {
            coef: vec![
                FiniteFieldElement::<ModEIP197>::from_value(6_i32), // constant term
                FiniteFieldElement::<ModEIP197>::from_value(-5_i32), // x term
                FiniteFieldElement::<ModEIP197>::from_value(1_i32), // x^2 term
            ],
        };
        assert_eq!(result, expected);
    }
}
