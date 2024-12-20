# Polynomials

### Definition: Polynomial

---

*A univariate polynomial over a commutative ring \\(R\\) with unity \\(1\\) is an expression of the form \\(f(x) = a_n x^{n} + a_{n-1} x^{n-1} + \cdots + a_1 x + a_0\\), where \\(x\\) is a variable and coefficients \\(a_0, \ldots, a_n\\) belong to \\(R\\). The set of all polynomials over \\(R\\) in the variable \\(x\\) is denoted as \\(R[x]\\).*

---

**Example:** In \\(\mathbb{Z}[x]\\), we have polynomials such as \\(2x^3 + x + 1\\), \\(x\\), and \\(1\\). In \\(\mathbb{R}[x]\\), we have polynomials like \\(\pi x^2 - \sqrt{2}x + e\\).

**Implementation:**

```rust
/// A struct representing a polynomial over a finite field.
#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial<F: Field> {
    /// Coefficients of the polynomial in increasing order of degree.
    pub coef: Vec<F>,
}
```

### Definition: Degree

---

*The degree of a non-zero polynomial \\(f(x) = a_n x^{n} + a_{n-1} x^{n-1} + \cdots + a_1 x + a_0\\), denoted as \\(\deg f\\), is the largest integer \\(n\\) such that \\(a_n \neq 0\\). The zero polynomial is defined to have degree \\(-1\\).*

---

**Example**

- \\(\deg(2x^3 + x + 1) = 3\\)
- \\(\deg(x) = 1\\)
- \\(\deg(1) = 0\\)
- \\(\deg(0) = -1\\)

**Implementation:**

```rust
impl<F: Field> Polynomial<F> {
    /// Removes trailing zeroes from a polynomial's coefficients.
    fn trim_trailing_zeros(poly: Vec<F>) -> Vec<F> {
        let mut trimmed = poly;
        while trimmed.last() == Some(&F::zero(None)) {
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
}
```

### Definition: Sum of Polynomials

---

*For polynomials \\(f(x) = \sum_{i=0}^n a_i x^i\\) and \\(g(x) = \sum_{i=0}^m b_i x^i\\), their sum is defined as: \\((f + g)(x) = \sum_{i=0}^{\max(n,m)} (a_i + b_i) x^i\\), where we set \\(a_i = 0\\) for \\(i > n\\) and \\(b_i = 0\\) for \\(i > m\\).*

---

**Example:** Let \\(f(x) = 2x^2 + 3x + 1\\) and \\(g(x) = x^3 - x + 4\\). Then, \\((f + g)(x) = x^3 + 2x^2 + 2x + 5\\)

**Implementation:**

```rust
impl<F: Field> Polynomial<F> {
    fn add_ref<'b>(&self, other: &'b Polynomial<F>) -> Polynomial<F> {
        let max_len = std::cmp::max(self.coef.len(), other.coef.len());
        let mut result = Vec::with_capacity(max_len);

        let zero = F::zero();

        for i in 0..max_len {
            let a = self.coef.get(i).unwrap_or(&zero);
            let b = other.coef.get(i).unwrap_or(&zero);
            result.push(a.add_ref(b));
        }
        Polynomial {
            coef: Self::trim_trailing_zeros(result),
        }
    }
}

impl<F: Field> Add for Polynomial<F> {
    type Output = Self;

    fn add(self, other: Self) -> Polynomial<F> {
        self.add_ref(&other)
    }
}
```

### Definition: Product of polynomials

---

*For polynomials \\(f(x) = \sum_{i=0}^n a_i x^i\\) and \\(g(x) = \sum_{j=0}^m b_j x^j\\), their product is defined as: \\((fg)(x) = \sum_{k=0}^{n+m} c_k x^k\\), where \\(c_k = \sum_{i+j=k} a_i b_j\\)*

---

**Example:** Let \\(f(x) = x + 1\\) and \\(g(x) = x^2 - 1\\). Then, \\((fg)(x) = x^3 + x^2 - x - 1\\)

**Implementation:**

```rust
impl<F: Field> Polynomial<F> {
    fn mul_ref<'b>(&self, other: &'b Polynomial<F>) -> Polynomial<F> {
        if self.is_zero() || other.is_zero() {
            return Polynomial::<F>::zero();
        }
        let mut result = vec![F::zero(); (self.degree() + other.degree() + 1) as usize];

        for (i, a) in self.coef.iter().enumerate() {
            for (j, b) in other.coef.iter().enumerate() {
                result[i + j] = result[i + j].add_ref(&a.mul_ref(b));
            }
        }
        Polynomial {
            coef: Polynomial::<F>::trim_trailing_zeros(result),
        }
    }
}

impl<F: Field> Mul<Polynomial<F>> for Polynomial<F> {
    type Output = Self;

    fn mul(self, other: Polynomial<F>) -> Polynomial<F> {
        self.mul_ref(&other)
    }
}
```

### Lemma 2.3.1

Let \\(K\\) be a field.

---

*Let \\(f, g \in K[x]\\) be non-zero polynomials. Then, \\(\deg(fg) = \deg f + \deg g\\).*

---

**Example:** Let \\(f(x) = x^2 + 1\\) and \\(g(x) = x^3 - x\\) in \\(\mathbb{R}[x]\\). Then, \\(\deg(fg) = \deg(x^5 - x^3 + x^2 + 1) = 5 = 2 + 3 = \deg f + \deg g\\)

### Theorem 2.3.2

We can also define division in the polynomial ring \\(K[x]\\).

---

*Let \\(f, g \in K[x]\\), with \\(g \neq 0\\). There exist unique polynomials \\(q, r \in K[x]\\) that satisfy \\(f = qg + r\\) and either \\(\deg r < \deg g\\) or \\(r = 0\\).*

---

**Proof**
    TBD

\\(q\\) is called the quotient of \\(f\\) divided by \\(g\\), and \\(r\\) is called the remainder; we write \\(r = f \bmod g\\).

**Example:** In \\(\mathbb{R}[x]\\), let \\(f(x) = x^3 + 2x^2 - x + 3\\) and \\(g(x) = x^2 + 1\\).  Then \\(f = qg + r\\) where \\(q(x) = x + 2\\) and \\(r(x) = -3x + 1\\).


**Implementation:**

```rust
impl<F: Field> Polynomial<F> {
    fn div_rem_ref<'b>(&self, other: &'b Polynomial<F>) -> (Polynomial<F>, Polynomial<F>) {
        if self.degree() < other.degree() {
            return (Polynomial::zero(), self.clone());
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
                remainder_coeffs[deg_diff + i] = remainder_coeffs[deg_diff + i]
                    .sub_ref(&(lead_term.mul_ref(&divisor_coeffs[i])));
            }
            remainder_coeffs = Self::trim_trailing_zeros(remainder_coeffs);
        }

        (
            Polynomial {
                coef: Self::trim_trailing_zeros(quotient),
            },
            Polynomial {
                coef: remainder_coeffs,
            },
        )
    }
}

impl<F: Field> Div for Polynomial<F> {
    type Output = Self;

    fn div(self, other: Polynomial<F>) -> Polynomial<F> {
        self.div_rem_ref(&other).0
    }
}
```

### Corollary

---

*Let \\(f \in K[x]\\) be a non-zero polynomial, and \\(a \in K\\) such that \\(f(a) = 0\\). Then, there exists a polynomial \\(q \in K[x]\\) such that \\(f(x) = (x - a)q(x)\\). In other words, \\((x - a)\\) is a factor of \\(f(x)\\).*

---

**Example:** Let \\(f(x) = x^2 + 1 \in (\mathbb{Z}/2\mathbb{Z})[x]\\). We have \\(f(1) = 1^2 + 1 = 0\\) in \\(\mathbb{Z}/2\mathbb{Z}\\), and indeed: \\(x^2 + 1 = (x - 1)^2 = x^2 - 2x + 1 = x^2 + 1\\) in \\((\mathbb{Z}/2\mathbb{Z})[x]\\)

### Theorem: Lagrange Interpolation

---

*A \\(n\\)-degre polynomial \\(P(x)\\) that goes through different \\(n + 1\\) points \\(\\{(x_1, y_1), (x_2, y_2), \cdots (x_{n + 1}, y_{n + 1})\\}\\) is uniquely represented as follows:*

\\[
    P(x) = \sum^{n+1}_{i=1} y_i \frac{f_i(x)}{f_i(x_i)}
\\]

*, where \\(f_i(x) = \prod_{k \neq i} (x - x_k)\\)*

---

**Proof** TBD

**Example:** The quadratic polynomial that goes through \\(\\{(1, 0), (2, 3), (3, 8)\\}\\) is as follows:

\begin{align*}
    P(x) = 0 \frac{(x - 2)(x - 3)}{(1 - 2) (1 - 3)} + 3 \frac{(x - 1)(x - 3)}{(2 - 1) (2 - 3)} + 8 \frac{(x - 1)(x - 2)}{(3 - 1) (3 - 2)} = x^{2} - 1
\end{align*}

Note that Lagrange interpolation finds the lowest degree of interpolating polynomial for the given vector.

**Implementation:**

```rust
impl<F: Field> Polynomial<F> {
    pub fn interpolate(x_values: &[F], y_values: &[F]) -> Polynomial<F> {
        let mut lagrange_polys = vec![];
        let numerators = Polynomial::from_monomials(x_values);

        for j in 0..x_values.len() {
            let mut denominator = F::one();
            for i in 0..x_values.len() {
                if i != j {
                    denominator = denominator * (x_values[j].sub_ref(&x_values[i]));
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
}
```

### Proposition: Homomorphisms of Lagrange Interpolation

*Let \\(L(v)\\) and \\(L(w)\\) be the polynomial resulting from Lagrange Interpolation on the output (\\(y\\)) vector \\(v\\) and \\(w\\) for the same inputs (\\(x\\)). Then, the following properties hold:*

- *Additivity: \\(L(v + w) = L(v) + L(w)\\) for any vectors \\(v\\) and \\(w\\)*
- *Scalar multiplication: \\(L(\gamma v) = \gamma L(v)\\) for any scalar \\(\gamma\\) and vector \\(v\\)*

**Proof**
    
Let \\(v = (v_1, \ldots, v_n)\\) and \\(w = (w_1, \ldots, w_n)\\) be vectors, and \\(x_1, \ldots, x_n\\) be the interpolation points.

\begin{align*}
    L(v + w) &= \sum_{i=1}^n (v_i + w_i) \prod_{j \neq i} \frac{x - x_j}{x_i - x_j} = \sum_{i=1}^n v_i \prod_{j \neq i} \frac{x - x_j}{x_i - x_j} + \sum_{i=1}^n w_i \prod_{j \neq i} \frac{x - x_j}{x_i - x_j} = L(v) + L(w) \\\\
    L(\gamma v) &= \sum_{i=1}^n (\gamma v_i) \prod_{j \neq i} \frac{x - x_j}{x_i - x_j} = \gamma \sum_{i=1}^n v_i \prod_{j \neq i} \frac{x - x_j}{x_i - x_j} = \gamma L(v)
\end{align*}
