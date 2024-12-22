# Arithmetization

The ultimate goal of Zero-Knowledge Proofs (ZKP) is to allow the prover to demonstrate their knowledge to the verifier without revealing any additional information. This knowledge typically involves solving complex problems, such as finding a secret input value that corresponds to a given public hash. ZKP protocols usually convert these statements into polynomial constraints. This process is often called **arithmetization**.

To make the protocol flexible, we need to encode this knowledge in a specific format, and one common approach is using Boolean circuits. It's well-known that problems in P (those solvable in polynomial time) and NP (those where the solution can be verified in polynomial time) can be represented as Boolean circuits. This means adopting Boolean circuits allows us to handle both P and NP problems.

However, Boolean circuits are often large and inefficient. Even a simple operation, like adding two 256-bit integers, can require hundreds of Boolean operators. In contrast, arithmetic circuits—essentially systems of equations involving addition, multiplication, and equality—offer a much more compact representation. Additionally, any Boolean circuit can be converted into an arithmetic circuit. For instance, the Boolean expression \\(z = x \land y\\) can be represented as \\(x(x-1) = 0\\), \\(y(y-1) = 0\\), and \\(z = xy\\) in an arithmetic circuit. Furthermore, as we'll see in this section, converting arithmetic circuits into polynomial constraints allows for much faster evaluation.

## Rank-1 Constraint System (R1CS)

There are many formats to represent arithmetic circuits, and one of the most popular ones is R1CS (Rank-1 Constraint System), which represents arithmetic circuits as \underline{a set of equality constraints, each involving only one multiplication}. In an arithmetic circuit, we call the concrete values assigned to the variables within the constraints witness. We first provide the formal definition of R1CS as follows:

### Definition: R1CS

An R1CS structure \\(\mathcal{S}\\) consists of:

- Size bounds \\(m, d, \ell \in \mathbb{N}\\) where \\(d > \ell\\)
- Three matrices \\(O, L, R \in \mathbb{F}^{m \times d}\\) with at most \\(\Omega(\max(m, d))\\) non-zero entries in total


An R1CS instance includes a public input \\(p \in \mathbb{F}^\ell\\), while an R1CS witness is a vector \\(w \in \mathbb{F}^{d - \ell - 1}\\).
A structure-instance tuple \\((S, p)\\) is satisfied by a witness \\(w\\) if:
\begin{equation}
(L \cdot a) \circ (R \cdot a) - O \cdot a = \mathbf{0}
\end{equation}
where \\(a = (1, w, p) \in \mathbb{F}^d\\), \\(\cdot\\) denotes matrix-vector multiplication, and \\(\circ\\) is the Hadamard product.

The intuitive interpretation of each matrix is as follows:


- \\(L\\): Encodes the left input of each gate
- \\(R\\): Encodes the right input of each gate
- \\(O\\): Encodes the output of each gate
- The leading 1 in the witness vector allows for constant terms


**Single Multiplication**

Let's consider a simple example where we want to prove \\(z = x \cdot y\\), with \\(z = 3690\\), \\(x = 82\\), and \\(y = 45\\).

- **Witness vector**: \\((1, z, x, y) = (1, 3690, 82, 45)\\)
- **Number of witnesses**: \\(m = 4\\)
- **Number of constraints**: \\(d = 1\\)


The R1CS constraint for \\(z = x \cdot y\\) is satisfied when:

\begin{align*}
&(\begin{bmatrix} 0 & 0 & 1 & 0 \end{bmatrix} \cdot a) \circ (\begin{bmatrix} 0 & 0 & 0 & 1 \end{bmatrix} \cdot a) - \begin{bmatrix} 0 & 1 & 0 & 0 \end{bmatrix} \cdot a \\\\
=&(\begin{bmatrix} 0 & 0 & 1 & 0 \end{bmatrix} \cdot \begin{bmatrix}
1 \\\\ 3690 \\\\ 82 \\\\ 45
\end{bmatrix}) \circ (\begin{bmatrix} 0 & 0 & 0 & 1 \end{bmatrix} \cdot \begin{bmatrix}
1 \\\\ 3690 \\\\ 82 \\\\ 45
\end{bmatrix}) - \begin{bmatrix} 0 & 1 & 0 & 0 \end{bmatrix} \cdot \begin{bmatrix}
1 \\\\ 3690 \\\\ 82 \\\\ 45
\end{bmatrix} \\\\
=& 82 \cdot 45 - 3690 \\\\
=& 3690 - 3690 \\\\
=& 0
\end{align*}

This example demonstrates how R1CS encodes a simple multiplication constraint:

- \\(L = \begin{bmatrix} 0 & 0 & 1 & 0 \end{bmatrix}\\) selects \\(x\\) (left input)
- \\(R = \begin{bmatrix} 0 & 0 & 0 & 1 \end{bmatrix}\\) selects \\(y\\) (right input)
- \\(O = \begin{bmatrix} 0 & 1 & 0 & 0 \end{bmatrix}\\) selects \\(z\\) (output)


**Multiple Constraints**

Let's examine a more complex example: \\(r = a \cdot b \cdot c \cdot d\\). Since R1CS requires that each constraint contain only one multiplication, we need to break this down into multiple constraints:

\begin{align*}
v_1 &= a \cdot b \\\\
v_2 &= c \cdot d \\\\
r &= v_1 \cdot v_2
\end{align*}

Note that alternative representations are possible, such as \\(v_1 = ab, v_2 = v_1c, r = v_2d\\). In this example, we use 7 variables \\((r, a, b, c, d, v_1, v_2)\\), so the dimension of the witness vector will be \\(m = 8\\) (including the constant 1). We have three constraints, so \\(n = 3\\).
To construct the matrices \\(L\\), \\(R\\), and \\(O\\), we can interpret the constraints as linear combinations:

\begin{align*}
v_1 &= (0 \cdot 1 + 0 \cdot r + 1 \cdot a + 0 \cdot b + 0 \cdot c + 0 \cdot d + 0 \cdot v_1 + 0 \cdot v_2) \cdot b \\\\
v_2 &= (0 \cdot 1 + 0 \cdot r + 0 \cdot a + 0 \cdot b + 1 \cdot c + 0 \cdot d + 0 \cdot v_1 + 0 \cdot v_2) \cdot d \\\\
r &= (0 \cdot 1 + 0 \cdot r + 0 \cdot a + 0 \cdot b + 0 \cdot c + 0 \cdot d + 1 \cdot v_1 + 0 \cdot v_2) \cdot v_2
\end{align*}

Thus, we can construct \\(L\\), \\(R\\), and \\(O\\) as follows:

\begin{equation*}
L = \begin{bmatrix}
0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\\\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\\\
0 & 0 & 0 & 0 & 0 & 0 & 1 & 0
\end{bmatrix}
\end{equation*}
\begin{equation*}
R = \begin{bmatrix}
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\\\
0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 \\\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1
\end{bmatrix}
\end{equation*}
\begin{equation*}
O = \begin{bmatrix}
0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \\\\
0 & 1 & 0 & 0 & 0 & 0 & 0 & 0
\end{bmatrix}
\end{equation*}

Where the columns in each matrix correspond to \\((1, r, a, b, c, d, v_1, v_2)\\).

**Addition with a Constant**

Let's examine the case \\(z = x \cdot y + 3\\). We can represent this as \\(-3 + z = x \cdot y\\). For the witness vector \\((1, z, x, y)\\), we have:

\begin{align*}
L &= \begin{bmatrix}
0 & 0 & 1 & 0
\end{bmatrix} \\\\
R &= \begin{bmatrix}
0 & 0 & 0 & 1
\end{bmatrix} \\\\
O &= \begin{bmatrix}
-3 & 1 & 0 & 0
\end{bmatrix}
\end{align*}

Note that the constant 3 appears in the \\(O\\) matrix with a negative sign, effectively moving it to the left side of the equation

**Multiplication with a Constant**

Now, let's consider \\(z = 3x^2 + y\\). The requirement of "one multiplication per constraint" doesn't apply to multiplication with a constant, as we can treat it as repeated addition.

\begin{align*}
L &= \begin{bmatrix}
0 & 0 & 3 & 0
\end{bmatrix} \\\\
R &= \begin{bmatrix}
0 & 0 & 1 & 0
\end{bmatrix} \\\\
O &= \begin{bmatrix}
0 & 1 & 0 & -1
\end{bmatrix}
\end{align*}

**Implementation:**

```rust
fn dot<F: Field>(a: &Vec<F>, b: &Vec<F>) -> F {
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
```

## Quadratic Arithmetic Program (QAP)

Recall that the prover aims to demonstrate knowledge of a witness \\(w\\) without revealing it. This is equivalent to knowing a vector \\(a\\) that satisfies \\((L \cdot a) \circ (R \cdot a) = O \cdot a\\), where \\(\circ\\) denotes the Hadamard (element-wise) product. However, evaluating this equivalence directly requires \\(\Omega(d)\\) operations, where \\(d\\) is the number of rows. To improve efficiency, we can convert this matrix comparison to a polynomial comparison, leveraging the Schwartz-Zippel Lemma, which allows us to check polynomial equality with \\(\Omega(1)\\) evaluations.

**Equivalence of Matrices**

Let's consider a simpler example to illustrate this concept. Suppose we want to test the equivalence \\(Av = Bu\\), where:

\begin{align*}
A = \begin{bmatrix}
2 & 5 \\\\
3 & 1 \\\\
\end{bmatrix},
B = \begin{bmatrix}
4 & 1 \\\\
2 & 3 \\\\
\end{bmatrix},
v = \begin{bmatrix}
3 \\\\ 1
\end{bmatrix},
u = \begin{bmatrix}
2 \\\\ 2
\end{bmatrix}
\end{align*}

The equivalence check can be represented as:

\begin{equation*}
\begin{bmatrix}
2 \\\\ 3
\end{bmatrix} \cdot 3 + \begin{bmatrix}
5 \\\\ 1
\end{bmatrix} \cdot 1 = \begin{bmatrix}
4 \\\\ 2
\end{bmatrix} \cdot 2 + \begin{bmatrix}
1 \\\\ 3
\end{bmatrix} \cdot 2
\end{equation*}

This matrix-vector equality check is equivalent to the following polynomial equality check:

\begin{equation*}
3 \cdot \lambda([(1, 2), (2, 3)]) + 1 \cdot \lambda([(1, 5), (2, 1)]) = 2 \cdot \lambda([(1, 4), (2, 2)]) + 2 \cdot \lambda([(1, 1), (2, 3)])
\end{equation*}

where \\(\lambda\\) denotes Lagrange Interpolation. In \\(\mathbb{F}_{11}\\) (field with 11 elements), we have:

\begin{align*}
\lambda([(1, 2), (2, 3)]) &= x + 1 \\\\
\lambda([(1, 5), (2, 1)]) &= 7x + 9 \\\\
\lambda([(1, 4), (2, 2)]) &= 9x + 6 \\\\
\lambda([(1, 1), (2, 3)]) &= 2x + 10
\end{align*}

The Schwartz-Zippel Lemma states that we need only one evaluation at a random point to check the equivalence of polynomials with high probability.

**Back to R1CS**

Let's thinks about how we can leverage the above method for the verification of R1CS. First, we can construct the interpolated polynomials for \\(L \cdot a\\), \\(R \cdot a\\), and \\(O \cdot a\\), denoted as \\(\ell(x)\\), \\(r(x)\\), and \\(o(x)\\), repectively, as follows:

\begin{align*}
\ell(x) = \sum^{d}_{i=1} a_i \ell_i(x) \quad \hbox{,where } \ell_i(x) := \lambda([(1, L_i,_1), (2, L_i,_2), \cdots,(m, L_i,_m)])
\end{align*}

\begin{align*}
r(x) = \sum^{d}_{i=1} a_i r_i(x) \quad \hbox{,where } r_i(x) := \lambda([(1, R_i,_1), (2, R_i,_2), \cdots,(m, R_i,_m)])
\end{align*}

\begin{align*}
o(x) = \sum^{d}_{i=1} a_i o_i(x) \quad \hbox{,where } o_i(x) := \lambda([(1, O_i,_1), (2, O_i,_2), \cdots,(m, O_i,_m)])
\end{align*}

However, the homomorphic property for multiplication doesn't hold for Lagrange Interpolation. While \\(\ell(x)\\), \\(r(x)\\), and \\(o(x)\\) are of degree at most \\(m-1\\), \\(\ell(x) \cdot r(x)\\) is of degree at most \\(2m-2\\). Thus, we don't have \\(\ell(x) \cdot r(x) = o(x)\\). 

To address this discrepancy, we introduce a degree \\(m\\) polynomial \\(t(x) = \prod_{i=1}^{m} (x - i)\\). Given the constituion of the interpolated equations, we have that \\(\forall{x} \in \\{1,\cdots, d\\}\\) \\(\ell(x) \cdot r(x) = o(x)\\). This implies the following:

\begin{equation}
\forall{x} \in \\{1,\cdots, m\\} \quad \ell(x) \cdot r(x) - o(x) = 0
\end{equation}

Thus, we can factorize \\(\ell(x) \cdot r(x) - o(x)\\) into the product of \\(t(x)\\) and an appripriate polynomial \\(h(x)\\) such that \\(\ell(x) \cdot r(x) - o(x) = t(x)h(x)\\).

Then, we can then rewrite the equation as:

\begin{equation}
\ell(x) \cdot r(x) = o(x) + h(x) \cdot t(x)
\end{equation}

This formulation allows us to maintain the desired polynomial relationships while accounting for the degree differences.

**Implementation:**

```rust
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
```
