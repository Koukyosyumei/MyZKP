# Gemini

[Bootle, Jonathan, et al. "Gemini: Elastic SNARKs for diverse environments." Annual International Conference on the Theory and Applications of Cryptographic Techniques. Cham: Springer International Publishing, 2022.](https://eprint.iacr.org/2022/420.pdf)

While many polynomial commitment schemes work for standard single-variable (univariate) polynomials, things get trickier when we deal with multilinear polynomials, which are polynomials with multiple variables where each variable has a maximum degree of one. The Gemini protocol offers an elegant way to handle these.

## Multilinear Polynomial

First, let's get a handle on the main object of interest. A multilinear polynomial with \\(m\\) variables (\\(f(X_0, X_1, \dots, X_{m- 1})\\)) can be written as a sum over all possible combinations of its variables. A general \\(m\\)-variate multilinear polynomial \\(f\\) is defined as:

\begin{align*}
    f(X_0, X_1, \dots, X_{m-1}) = \sum_{i = 0}^{2^m-1} c_i \prod_{j=0}^{m-1} X_j^{b_j(i)}
\end{align*}

Here, \\(c_i\\) is the \\(i\\)-th coefficient, and \\((b_{m-1}(i), \dots, b_0(i))\\) is the binary representation of the number \\(i\\). THis formula just means we sum up terms where each variable \\(X_{j}\\) is either present or absent.

For example, a three-variable (\\(m = 3\\)) multilinear polynomial looks like this:

\begin{align*}
    f(X_0, X_1, X_2) = c_{000} + c_{100} X_0 + c_{010} X_1 + c_{110} X_0 X_1 + c_{001} X_2 + c_{101} X_0 X_2 + c_{011} X_1 X_2 + c_{111} X_0 X_1 X_2
\end{align*}

The subscript of each coefficient \\(c\\) tells you which variables are in that term. For instance, \\(c_{110}\\) is the coefficient of the \\(X_0 X_1\\) term (reading the bits from right to left as \\(i_0\\),\\(i_1\\),\\(i_2\\)).

## Problem

Our goal is to create a commitment scheme for these polynomials. Specifically, a Prover wants to commit to a multilinear polynomial \\(f\\) and then convince a Verifier that for a specific, publicly known point \\(\rho = (\rho_0, \rho_1, \dots, \rho_{m-1})\\) and a public value \\(u\\), the following statement is true:

\begin{align*}
    f(\rho) = u
\end{align*}

The Prover must prove this without revealing the entire polynomial \\(f\\).

## Key idea 

Here's the core idea of Gemini: **transform the multivariable problem into a series of single-variable problems**, which are much easier to handle.

We first take the \\(2^{m}\\) coefficients of our multilinear polynomial \\(f\\) and use them to define a new **univariate** polynomial \\(g\\) of degree \\(2^{m} - 1\\). 

\begin{align*}
    g(X) = \sum_{i=0}^{2^{m} - 1} c_{i} X^{i}
\end{align*}

Using our three-variable example, the corresponding univariate polynomial \\(g(X)\\) would be:

\begin{align*}
g(X) = c_{000} + c_{100} X + c_{010} X^2 + c_{110} X^3 + c_{001} X^4 + c_{101} X^5 + c_{011} X^6 + c_{111} X^7
\end{align*}

### Split and Fold

Then, we are going to recursively fold this univariate polynomialan univariate polynomial with the **folding** trick, which repeatedly halves the degree by separating even and odd exponents:

For any univariate polynomial \\(g(X)\\), we can write:

\begin{align*}
    g(X) = g_E(X) + X \cdot g_O(X)
\end{align*}, where

\begin{align*}
    g_E(X) = \frac{g(X) + g(-X)}{2}, \quad \quad g_O(X) = \frac{g(X) - g(-X)}{2X}
\end{align*}

Here, \\(g_E\\) contains the even exponents of \\(g\\) and \\(g_O\\) contains the shifted odd exponents.

\begin{align*}
    g_E(X) = \sum_{i=0}^{\frac{m+1}{2} - 1} c_{2i} X^{2i}, \quad \quad g_O(X) = \sum_{i=0}^{\frac{m+1}{2} - 1} c_{2i + 1} X^{2i}
\end{align*}

Now treat \\(Y := X^2\\). Because the even and odd parts only use powers \\(X^{2k}\\), both \\(g_E\\) and \\(g_O\\) can be seen as polynomials in \\(Y\\) with half the degree.

We then define a sequence of folded polynomials as follows:

\begin{align*}
    g^{(0)}(X) &:= g(X) \\\\
    g^{(1)}(X) &:= g_E^{(0)}(X) + \rho_0 g_O^{(0)}(X) \\\\
    g^{(2)}(X^2) &:= g_E^{(1)}(X^2) + \rho_1 g_O^{(1)}(X^2) \\\\
        &\quad\vdots \\\\
    g^{(m)}(X^{2^{m-1}}) &:= g_E^{(m-1)}(X^{2^{m-1}}) + \rho_{m-1} g_O^{(m-1)}(X^{2^{m-1}}) = (\text{a constant}) \\\\
\end{align*}

Intuitively, each folding step substitutes the known scalar \\(\rho_{i-1}\\) for one of the original variables, and reduces the number of remaining variables by one. After \\(m\\) folds, we obtain a constant equal to \\(f(\rho)\\), i.e., \\(g^{(m)}(\cdot) = u\\).

```rust
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
```

For example, the above three-variable example yeilds the following sequence, where \\(Y_i\\) denotes \\(X^{2^{i}}\\):

\begin{align*}
g^{(0)}(Y_0) &= c_{000} + c_{100} Y_0 + c_{010} Y_0^2 + c_{110} Y_0^3 + c_{001} Y_0^4 + c_{101} Y_0^5 + c_{011} Y_0^6 + c_{111} Y_0^7 \\\\
             &= c_{000} + c_{010} Y_0^2 + c_{001} Y_0^4 + c_{011} Y_0^6 + Y_0 (c_{100} + c_{110} Y_0^2 + c_{101} Y_0^4 + c_{111} Y_0^6) \\\\
g^{(1)}(Y_1) &= g_{E}^{(0)}(Y_1) + \rho_0 g_{O}^{(0)}(Y_1) \\\\ 
             &= c_{000} + c_{010} Y_1 + c_{001} Y_1^2 + c_{011} Y_1^3 + \rho_0 (c_{100} + c_{110} Y_1 + c_{101} Y_1^2 + c_{111} Y_1^3) \\\\
             &= (c_{000} + \rho_0  c_{100}) + (c_{010} + \rho_0  c_{110}) Y_1 + (c_{001} + \rho_0  c_{101}) Y_1^2 + (c_{011} + \rho_0  c_{111}) Y_1^3 \\\\
             &= (c_{000} + \rho_0  c_{100}) + (c_{001} + \rho_0  c_{101}) Y_1^2 + Y_1 ((c_{010} + \rho_0  c_{110}) + (c_{011} + \rho_0  c_{111}) Y_1^2) \\\\
g^{(2)}(Y_2) &= g_{E}^{(1)}(Y_2) + \rho_1 g_{O}^{(1)}(Y_2) \\\\ 
             &= (c_{000} + \rho_0  c_{100}) + (c_{001} + \rho_0  c_{101}) Y_2 + \rho_1 ((c_{010} + \rho_0  c_{110}) + (c_{011} + \rho_0  c_{111}) Y_2) \\\\
             &= (c_{000} + \rho_0  c_{100} + \rho_1 c_{010} + \rho_0 \rho_1 c_{110}) + Y_2 (c_{001} + \rho_0 c_{101} + \rho_1 c_{011} + \rho_0 \rho_1 c_{111}) \\\\
g^{(3)}(Y_3) &= g_{E}^{(1)}(Y_3) + \rho_2 g_{O}^{(1)}(Y_3) \\\\
             &= (c_{000} + \rho_0  c_{100} + \rho_1 c_{010} + \rho_0 \rho_1 c_{110}) + \rho_2 (c_{001} + \rho_0 c_{101} + \rho_1 c_{011} + \rho_0 \rho_1 c_{111}) \\\\
             &= c_{000} + \rho_0  c_{100} + \rho_1 c_{010} + \rho_0 \rho_1 c_{110} + \rho_2 c_{001} + \rho_0 \rho_2 c_{101} + \rho_1 \rho_2 c_{011} + \rho_0 \rho_1 \rho_2 c_{111} \\\\
             &= u
\end{align*}

### Relationship used for verification

The recurrence used between consecutive folded polynomials is:

\begin{align*}
    g^{(i)}(X^2) &= g_{E}^{(i-1)}(X) + \rho_{i-1} g_{O}^{(i-1)}(X) \\\\ 
               &=\frac{g^{(i-1)}(X) + g^{(i-1)}(-X)}{2} + \rho_{i-1} \frac{g^{(i-1)}(X) - g^{(i-1)}(-X)}{2X}
\end{align*}

By the Schwartz–Zippel lemma, we can check this relationship by evaluating the polynomials at a random point, with only negligible probability of successful cheating.

## The protocol

At a high level, Gemini works as follows:

- **Pubic Inputs**: the evaluation points \\(\rho\\) and the value \\(u\\)
- **Goal**: prover convinces verifier that \\(f(\rho) = u\\) (equivalently \\(g^{(m)}(\cdot) = u\\))

1. The prover commits all \\(g^{(i)}\\)
2. Verifier samples a random value \\(\beta\\)
3. For each round \\(u = 1, \dots, m\\), the verifier asks the prover to open:
   - \\(g^{i-1}\\) at \\(\beta\\),
   - \\(g^{i-1}\\) at \\(-\beta\\),
   - \\(g^{i}\\) at \\(\beta^2\\)
4. The verifier checks, for each \\(i\\), that the three opened values satisfy the recurrence:

\begin{align*}
    g^{(i)}(\beta^2) = \frac{g^{(i-1)}(\beta) + g^{(i-1)}(-\beta)}{2} + \rho_{i-1} \frac{g^{(i-1)}(\beta) - g^{(i-1)}(-\beta)}{2\beta}
\end{align*}

5. Finally, the verifier checks that the last folded value \\(g^{(m)}(\beta^2)\\) is equal to the public target \\(u\\).

For example, we can use KZG commitment scheme to commit and open the above polynomials. We can leverage the batch proof of KZG to evaluate a polynomial at \\(\beta\\), \\(-\beta\\), and \\(\beta^2\\) at once, and degree-bound test to check the degree of each polynomial.

```rust
pub type CommitmentGemini = Vec<CommitmentKZG>;

pub struct ProofGemini {
    es: Vec<ProofKZG>,
    es_neg: Vec<ProofKZG>,
    es_hat: Vec<ProofKZG>,
    degree_proofs: Vec<ProofDegreeBound>,
}

pub fn commit_gemini(polys: &[Polynomial<FqOrder>], pk: &PublicKeyKZG) -> CommitmentGemini {
    polys.iter().map(|p| commit_kzg(p, pk)).collect()
}

pub fn open_gemini(
    polys: &[Polynomial<FqOrder>],
    beta: &FqOrder,
    pk: &PublicKeyKZG,
) -> ProofGemini {
    let num_polys = polys.len();
    ProofGemini {
        es: polys
            .iter()
            .take(polys.len() - 1)
            .map(|p| open_kzg(p, beta, pk))
            .collect(),
        es_neg: polys
            .iter()
            .take(polys.len() - 1)
            .map(|p| open_kzg(p, &(FqOrder::zero() - beta).sanitize(), pk))
            .collect(),
        es_hat: polys
            .iter()
            .skip(1)
            .take(polys.len() - 2)
            .map(|p| open_kzg(p, &beta.pow(2), pk))
            .collect(),
        degree_proofs: polys
            .iter()
            .enumerate()
            .map(|(i, p)| prove_degree_bound(p, pk, 2_usize.pow((num_polys - i - 1) as u32)))
            .collect(),
    }
}

pub fn verify_gemini(
    rhos: &Vec<FqOrder>,
    mu: &FqOrder,
    beta: &FqOrder,
    commitment: &CommitmentGemini,
    proof: &ProofGemini,
    pk: &PublicKeyKZG,
) -> bool {
    let log2_n = rhos.len();
    if log2_n != commitment.len() - 1 {
        return false;
    }

    if !commitment
        .iter()
        .zip(proof.degree_proofs.iter())
        .enumerate()
        .all(|(i, (c, p))| verify_degree_bound(c, p, pk, 2_usize.pow((log2_n - i) as u32)))
    {
        return false;
    }

    if !commitment
        .iter()
        .take(commitment.len() - 1)
        .zip(proof.es.iter())
        .all(|(c, p)| verify_kzg(beta, c, p, pk))
    {
        return false;
    }

    if !commitment
        .iter()
        .take(commitment.len() - 1)
        .zip(proof.es_neg.iter())
        .all(|(c, p)| verify_kzg(&(FqOrder::zero() - beta).sanitize(), c, p, pk))
    {
        return false;
    }

    if !commitment
        .iter()
        .skip(1)
        .take(commitment.len() - 2)
        .zip(proof.es_hat.iter())
        .all(|(c, p)| verify_kzg(&beta.pow(2), c, p, pk))
    {
        return false;
    }
    let es = proof.es.iter().map(|p| p.y.clone()).collect::<Vec<_>>();
    let es_neg = proof.es_neg.iter().map(|p| p.y.clone()).collect::<Vec<_>>();
    let mut es_hat = proof.es_hat.iter().map(|p| p.y.clone()).collect::<Vec<_>>();
    es_hat.push(mu.clone());

    let two = FqOrder::from_value(2);
    (0..log2_n).all(|j| {
        two.mul_ref(beta).mul_ref(&es_hat[j])
            == beta.mul_ref(&es[j].add_ref(&es_neg[j]))
                + rhos[j].mul_ref(&es[j].sub_ref(&es_neg[j]))
    })
}
```

### Example

Take

\begin{align*}
    f(X_0, X_1, X_2) = 1 + 2 X_0 + 3 X_1 + 4 X_0 X_1 + 5 X_2 + 6 X_0 X_2 + 7 X_1 X_2 + 8 X_0 X_1 X_2
\end{align*}

so when packed into the univariate \\(g(X) = 1 + 2 X + 3 X^2 + 4 X^3 + 5 X^4 + 6 X^5 + 7 X^6 + 8 X^7\\).

Let \\(rho = (1, 2, 3)\\). Then, \\f(\rho) = 140\\, which will be out public target \\(u\\) and yeild the following folded polynomials:

\begin{align*}
    g^{(0)}(Y_0) &= 1 + 3 Y_0^2 + 5 Y_0^4 + 7 Y_0^6 + Y_0 (2 + 4 Y_0^2 + 6 Y_0^4 + 8 Y_0^6)  \\\\
    g^{(1)}(Y_1) &= 1 + 3 Y_1 + 5 Y_1^2 + 7 Y_1^3 + (2 + 4 Y_1 + 6 Y_1^2 + 8 Y_1^3) \\\\
                 &= 3 + 7 Y_1 + 11 Y_1^2 + 15 Y_1^3\\\\
                 &= 3 + 11 Y_1^2 + Y_1 (7 + 15 Y_1^2) \\\\
    g^{(2)}(Y_2) &= 3 + 11 Y_2 + 2 (7 + 15 Y_2)\\\\
                 &= 17 + Y_2 (41) \\\\
    g^{(3)}(Y_3) &= 17 + 3 \cdot 41 = 140
\end{align*}

Choose \\(\beta = 2\\). Then, we can observe that the fold sequence evaluates numerically as:

\begin{align*}
    g^{(0)}(2) &= 1793, \quad g^{(0)}(-2) = -711, \quad g^{(1)}(4) = 1167, \quad \frac{g^{(0)}(2) + g^{(0)}(-2)}{2} + 1 \cdot \frac{g^{(0)}(2) - g^{(0)}(-2)}{2 * 2} = 1167 \\\\
    g^{(1)}(2) &= 181, \quad g^{(1)}(-2) = -87, \quad g^{(2)}(4) = 181, \quad \frac{g^{(1)}(2) + g^{(1)}(-2)}{2} + 2 \cdot \frac{g^{(1)}(2) - g^{(1)}(-2)}{2 * 2} = 181 \\\\
    g^{(2)}(2) &= 99, \quad g^{(2)}(-2) = -65, \quad g^{(3)}(4) = u = 140, \quad \frac{g^{(2)}(2) + g^{(2)}(-2)}{2} + 3 \cdot \frac{g^{(2)}(2) - g^{(2)}(-2)}{2 * 2} = 140
\end{align*}
 