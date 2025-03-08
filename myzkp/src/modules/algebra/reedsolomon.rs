use lazy_static::lazy_static;
use num_bigint::BigInt;
use num_traits::{One, Zero};
use paste::paste;
use serde::{Deserialize, Serialize};
use std::iter;
use std::str::FromStr;

use crate::define_extension_field;
use crate::define_myzkp_modulus_type;

use crate::modules::algebra::efield::ExtendedFieldElement;
use crate::modules::algebra::efield::IrreduciblePoly;
use crate::modules::algebra::field::Field;
use crate::modules::algebra::field::FiniteFieldElement;
use crate::modules::algebra::field::ModulusValue;
use crate::modules::algebra::polynomial::Polynomial;
use crate::modules::algebra::ring::Ring;

pub struct ReedSolomon<F: Field> {
    n: usize, // codeword length
    d: usize, // validity length
    k: usize, // message length
    g: F,     // generator
}

impl<F: Field> ReedSolomon<F> {
    pub fn new(n: usize, k: usize, g: F) -> Self {
        assert!(n >= k, "n must be at least k");
        let d = n - k;
        ReedSolomon { n, d, k, g }
    }

    pub fn generator_polynomial(&self) -> Polynomial<F> {
        let eval_points = self.evaluation_points(self.d);
        Polynomial::<F>::from_monomials(&eval_points)
    }

    /// Compute evaluation points: [g^0, g^1, ..., g^{el - 1}]
    fn evaluation_points(&self, el: usize) -> Vec<F> {
        let mut points = Vec::with_capacity(el);
        let one = F::one();
        for i in 0..el {
            let pt = self.g.pow(i);
            points.push(pt);
        }
        points
    }

    /// Encode a message (as a slice of field elements representing polynomial coefficients,
    /// from constant term up) into a codeword of length n.
    ///
    /// The message polynomial must have degree less than k.
    pub fn encode(&self, message: &[F]) -> Vec<F> {
        assert!(
            message.len() <= self.k,
            "Received message must be less than or equal to k (k={}, get {})",
            self.k,
            message.len()
        );

        let m_poly = Polynomial {
            coef: message.to_vec(),
        };

        // let shift_length = self.d - message.len();
        let mut m_shifted_coef = vec![F::from_value(0); self.d];
        m_shifted_coef.extend(m_poly.coef);

        let m_shifted = Polynomial {
            coef: m_shifted_coef,
        };

        let g = self.generator_polynomial();
        let remainder = &m_shifted % &g;

        let mut codeword = m_shifted - remainder;
        codeword.coef
    }

    /// Compute syndromes S₁, S₂, …, S₂ₜ from the received codeword.
    /// Here t = (n - k) / 2 and d = n - k.
    pub fn compute_syndromes(&self, received: &[F]) -> Vec<F> {
        let eval_points = self.evaluation_points(self.n);
        let mut syndromes = Vec::with_capacity(self.d);
        for j in 0..self.d {
            let mut s = F::zero();
            for (i, r_i) in received.iter().enumerate() {
                // (α^i)^j = α^(i * j)
                s = s + r_i.clone() * eval_points[i].clone().pow(j as i32);
            }
            syndromes.push(s);
        }
        syndromes
    }

    /// Berlekamp–Massey algorithm: given syndromes, compute the error locator polynomial σ(x).
    /// The returned polynomial is of the form σ(x) = 1 + σ₁x + … + σₜxᵗ.
    fn berlekamp_massey(&self, syndromes: &[F]) -> Polynomial<F> {
        let n = syndromes.len();
        let mut sigma = Polynomial {
            coef: vec![F::one()],
        };
        let mut B = Polynomial {
            coef: vec![F::one()],
        };
        let mut L = 0;
        let mut m = 1;
        let mut b = F::one();

        for n_iter in 0..n {
            // Compute discrepancy d = S[n_iter] + Σ₍ᵢ₌₁₎^L σᵢ S[n_iter - i]
            let mut d = syndromes[n_iter].clone();
            for i in 1..=L {
                if i < sigma.coef.len() {
                    d = d + sigma.coef[i].clone() * syndromes[n_iter - i].clone();
                }
            }
            if d == F::zero() {
                m += 1;
            } else {
                let T = sigma.clone();
                // factor = d / b
                let factor = d.clone() / b.clone();
                // Multiply B(x) by x^m: create x^m with coefficient 1 in degree m.
                let x_m = Polynomial {
                    coef: iter::repeat(F::zero())
                        .take(m)
                        .chain(iter::once(F::one()))
                        .collect(),
                };
                let x_m_B = x_m * B.clone();
                // sigma = sigma - factor * (x^m * B(x))
                sigma = sigma - x_m_B * factor.clone();
                if 2 * L <= n_iter {
                    L = n_iter + 1 - L;
                    B = T;
                    b = d;
                    m = 1;
                } else {
                    m += 1;
                }
            }
        }
        sigma
    }

    /// Compute the formal derivative of a polynomial.
    fn polynomial_derivative(poly: &Polynomial<F>) -> Polynomial<F> {
        let deriv_coeffs = poly
            .coef
            .iter()
            .enumerate()
            .skip(1) // skip the constant term
            .map(|(i, coef)| coef.clone() * F::from_value(i as i32))
            .collect();
        Polynomial { coef: deriv_coeffs }
    }

    /// Compute the error evaluator polynomial Ω(x) = [σ(x) * S(x)] mod x^(2t),
    /// where S(x) is the syndrome polynomial.
    fn compute_error_evaluator(&self, syndromes: &[F], sigma: &Polynomial<F>) -> Polynomial<F> {
        let t = (self.n - self.k) / 2;
        // Build syndrome polynomial S(x) = S₁ + S₂x + … + S₂ₜ x^(2t-1).
        let S_poly = Polynomial {
            coef: syndromes.to_vec(),
        };
        let prod = sigma.clone() * S_poly;
        // Truncate the product to degree < 2t.
        let mut truncated = prod;
        if truncated.coef.len() > 2 * t {
            truncated.coef.truncate(2 * t);
        }
        truncated
    }

    /// Find error locations by searching for roots of σ(x).
    ///
    /// In our RS code the error locator polynomial is defined so that if an error occurred at
    /// position i (with evaluation point α^i), then x = (α^i)⁻¹ is a root.
    /// We return the indices where an error is detected.
    fn find_error_locations(&self, sigma: &Polynomial<F>) -> Vec<usize> {
        let mut error_positions = Vec::new();
        let eval_points = self.evaluation_points(self.n);
        for (i, pt) in eval_points.iter().enumerate() {
            // Compute x = (α^i)⁻¹.
            let pt_inv = pt.clone().inverse();
            if sigma.eval(&pt_inv) == F::zero() {
                error_positions.push(i);
            }
        }
        error_positions
    }

    /// Correct errors in a received codeword (with unknown errors) using syndrome
    /// computation, the Berlekamp–Massey algorithm, and Forney’s formula.
    ///
    /// Returns the corrected codeword if error correction succeeds, or None if there are too many errors.
    pub fn decode_errors(&self, received: &[F]) -> Option<Vec<F>> {
        assert!(
            received.len() == self.n,
            "Received codeword must have length n"
        );
        let syndromes = self.compute_syndromes(received);
        // If all syndromes are zero, no error occurred.
        if syndromes.iter().all(|s| *s == F::zero()) {
            return Some(received.to_vec());
        }
        // Compute the error locator polynomial σ(x).
        let sigma = self.berlekamp_massey(&syndromes);
        // Find error positions.
        let error_positions = self.find_error_locations(&sigma);
        // The expected number of errors is the degree of σ(x).
        let num_errors = sigma.degree();
        if error_positions.len() != num_errors.try_into().unwrap() {
            // The number of roots found does not match the degree (decoding failure).
            return None;
        }
        // Compute error evaluator polynomial Ω(x).
        let omega = self.compute_error_evaluator(&syndromes, &sigma);
        // Compute the derivative σ'(x).
        let sigma_deriv = Self::polynomial_derivative(&sigma);

        // Correct errors using Forney’s formula.
        let mut corrected = received.to_vec();
        let eval_points = self.evaluation_points(self.n);
        for &pos in error_positions.iter() {
            let Xi = eval_points[pos].clone();
            let Xi_inv = Xi.clone().inverse();
            let omega_val = omega.eval(&Xi_inv);
            let sigma_deriv_val = sigma_deriv.eval(&Xi_inv);
            if sigma_deriv_val == F::zero() {
                // Cannot correct if derivative evaluates to zero.
                return None;
            }
            // Forney’s formula: error magnitude = - (Xi * Ω(Xi⁻¹)) / σ'(Xi⁻¹)
            let error_mag = -(Xi * omega_val) / sigma_deriv_val;
            // Correct the error (assuming r = c + e, we subtract the error magnitude).
            corrected[pos] = corrected[pos].clone() - error_mag;
        }
        // Optionally, one could recompute syndromes on the corrected word.
        Some(corrected)
    }
}

define_myzkp_modulus_type!(Mod2, "2");
define_extension_field!(
    Ip0x11D,
    FiniteFieldElement<Mod2>,
    Polynomial {
        coef: vec![
            FiniteFieldElement::<Mod2>::from_value(1),
            FiniteFieldElement::<Mod2>::from_value(0),
            FiniteFieldElement::<Mod2>::from_value(1),
            FiniteFieldElement::<Mod2>::from_value(1),
            FiniteFieldElement::<Mod2>::from_value(1),
            FiniteFieldElement::<Mod2>::from_value(0),
            FiniteFieldElement::<Mod2>::from_value(0),
            FiniteFieldElement::<Mod2>::from_value(0),
            FiniteFieldElement::<Mod2>::from_value(1),
        ],
    }
);

type GF2to8 = ExtendedFieldElement<Mod2, Ip0x11D>;
impl GF2to8 {
    pub fn from_u8(n: u8) -> Self {
        let mut coef = Vec::with_capacity(8);
        for i in 0..8 {
            let bit = (n >> i) & 1;
            coef.push(FiniteFieldElement::<Mod2>::from_value(bit as i32));
        }
        GF2to8::new(Polynomial { coef })
    }
}

pub struct ReedSolomon2D<F: Field> {
    col_coder: ReedSolomon<F>,
    row_coder: ReedSolomon<F>,
}

impl<F: Field> ReedSolomon2D<F> {
    pub fn new(col_codeword_len: usize, row_codeword_len: usize, message_len: usize, g: F) -> Self {
        let col_coder = ReedSolomon::new(col_codeword_len, message_len, g.clone());
        let row_coder = ReedSolomon::new(row_codeword_len, message_len, g.clone());
        ReedSolomon2D {
            col_coder,
            row_coder,
        }
    }

    pub fn encode(&self, data: &[F]) -> Vec<Vec<F>> {
        let mut matrix = self.organize_data_matrix(data);

        // First dimension encoding (rows)
        // - message_len x row_codeword_len
        let encoded_row = matrix
            .iter()
            .map(|row| self.row_coder.encode(row))
            .collect::<Vec<_>>();

        // Second dimension encoding (columns)
        // - row_codeword_len x message_len
        let transposed_encoded_row = self.transpose_matrix(&encoded_row);
        // - row_codeword_len x col_codeword_len
        let encoded_cols = transposed_encoded_row
            .iter()
            .map(|col| self.col_coder.encode(col))
            .collect::<Vec<_>>();

        // col_codeword_len x row_codeword_len
        self.transpose_matrix(&encoded_cols)
    }

    fn organize_data_matrix(&self, data: &[F]) -> Vec<Vec<F>> {
        let size = (data.len() as f64).sqrt().ceil() as usize;
        let mut matrix = vec![vec![F::zero(); size]; size];

        data.iter().enumerate().for_each(|(i, val)| {
            let row = i / size;
            let col = i % size;
            matrix[row][col] = val.clone();
        });

        matrix
    }

    fn transpose_matrix(&self, matrix: &[Vec<F>]) -> Vec<Vec<F>> {
        let row_size = matrix[0].len();
        let col_size = matrix.len();
        let mut transposed = vec![vec![F::zero(); col_size]; row_size];
        for (row_idx, row) in matrix.iter().enumerate() {
            for (col_idx, val) in row.iter().enumerate() {
                transposed[col_idx][row_idx] = val.clone();
            }
        }
        transposed
    }
}

pub fn create_rs(codeword_len: usize, message_len: usize) -> ReedSolomon<GF2to8> {
    let coef = vec![
        FiniteFieldElement::<Mod2>::zero(),
        FiniteFieldElement::<Mod2>::one(),
    ];
    let generator = GF2to8::new(Polynomial { coef });
    ReedSolomon::new(codeword_len, message_len, generator)
}

pub fn create_rs2d(
    col_codeword_len: usize,
    row_codeword_len: usize,
    message_len: usize,
) -> ReedSolomon2D<GF2to8> {
    let coef = vec![
        FiniteFieldElement::<Mod2>::zero(),
        FiniteFieldElement::<Mod2>::one(),
    ];
    let generator = GF2to8::new(Polynomial { coef });
    ReedSolomon2D::new(col_codeword_len, row_codeword_len, message_len, generator)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_syndrome_computation_no_error() {
        let rs = create_rs(7, 3);
        let message = vec![GF2to8::from_u8(9), GF2to8::from_u8(1), GF2to8::from_u8(7)];
        let codeword = rs.encode(&message);
        let syndromes = rs.compute_syndromes(&codeword);
        for (j, s) in syndromes.iter().enumerate() {
            assert_eq!(
                *s,
                GF2to8::from_u8(0),
                "Syndrome S_{} should be zero for an error-free codeword",
                j + 1
            );
        }
    }

    /// Test that when no errors occur, decode_errors returns the original codeword.
    #[test]
    fn test_no_errors() {
        let rs = create_rs(7, 3);
        let message = vec![GF2to8::from_u8(9), GF2to8::from_u8(1), GF2to8::from_u8(7)];
        let codeword = rs.encode(&message);
        println!("aaa {}", codeword.len());
        let decoded = rs
            .decode_errors(&codeword)
            .expect("Decoding should succeed with no errors");
        assert_eq!(
            codeword, decoded,
            "Decoding an error-free codeword should yield the original codeword"
        );
    }

    #[test]
    fn test_single_error_correction() {
        let rs = create_rs(7, 3);
        let message = vec![GF2to8::from_u8(1), GF2to8::from_u8(2), GF2to8::from_u8(3)];
        let mut codeword = rs.encode(&message);
        // Introduce a single error at position 3: add an offset to the symbol.
        codeword[3] = codeword[3].clone() + GF2to8::from_u8(10);
        let corrected = rs
            .decode_errors(&codeword)
            .expect("Decoding should correct a single error");
        let expected = rs.encode(&message);
        assert_eq!(
            corrected, expected,
            "After correcting a single error, the codeword should match the original encoding"
        );
    }

    #[test]
    fn test_multiple_error_correction() {
        let rs = create_rs(7, 3);
        let message = vec![GF2to8::from_u8(4), GF2to8::from_u8(8), GF2to8::from_u8(15)];
        let mut codeword = rs.encode(&message);
        // Introduce errors in two positions.
        codeword[1] = codeword[1].clone() + GF2to8::from_u8(7);
        codeword[5] = codeword[5].clone() + GF2to8::from_u8(12);
        let corrected = rs
            .decode_errors(&codeword)
            .expect("Decoding should correct two errors");
        let expected = rs.encode(&message);
        assert_eq!(
            corrected, expected,
            "After correcting two errors, the codeword should match the original encoding"
        );
    }

    #[test]
    fn test_too_many_errors() {
        let rs = create_rs(7, 3);
        let message = vec![GF2to8::from_u8(2), GF2to8::from_u8(4), GF2to8::from_u8(6)];
        let mut codeword = rs.encode(&message);
        // For n = 7 and k = 3, the decoder can correct up to 2 errors.
        // Introduce three errors (exceeding the correction capability).
        codeword[0] = codeword[0].clone() + GF2to8::from_u8(3);
        codeword[3] = codeword[3].clone() + GF2to8::from_u8(5);
        codeword[6] = codeword[6].clone() + GF2to8::from_u8(7);
        let decoded = rs.decode_errors(&codeword);
        assert!(
            decoded.is_none(),
            "Decoding should fail when more errors occur than the code can correct"
        );
    }

    #[test]
    fn test_no_errors_2d() {
        let rs2d = create_rs2d(7, 7, 3);
        let message = vec![GF2to8::from_u8(2), GF2to8::from_u8(4), GF2to8::from_u8(6)];
        let codeword = rs2d.encode(&message);
    }
}
