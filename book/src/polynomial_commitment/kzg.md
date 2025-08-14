# KZG

Let \\(G_1\\), \\(G_2\\), and \\(G_T\\) be groups such that there exists a bilinear mapping \\(e: G_1 \times G_2 \to G_T\\). Let \\(g_1 \in G_1\\) and \\(g_2 \in G_2\\) be fixed generators.

We want a short commitment to a polynomial over \\(F\\) denoted as \\(f(X) = \sum_{i=0}^{d} c_{i} x^{i} \in F_{q}[X]\\) with a public maximum degree bound \\(D\\) such that \\(d \leq D\\).

## Setup

A trusted party samples a secret \\(\alpha \xleftarrow{\text{random}} F_{q}\\), and publishes a Structured Reference String (SRS):

\begin{align*}
    &\mathrm{PK} := (\\{g_1^{(\alpha^{i})}\\}^{D}_{i=0}, g_2, g^{\alpha}_2)
\end{align*}

```rust
pub fn setup_kzg(g1: &G1Point, g2: &G2Point, n: usize) -> PublicKeyKZG {
    let alpha = FqOrder::random_element(&[]);

    let mut alpha_1 = Vec::with_capacity(n);
    let mut s_power = FqOrder::one();
    for _ in 0..1 + n {
        alpha_1.push(g1.mul_ref(s_power.clone().get_value()));
        s_power = s_power * alpha.clone();
    }

    let alpha_2 = vec![g2.clone(), g2.mul_ref(alpha.get_value())];

    PublicKeyKZG { alpha_1, alpha_2 }
}
```

## Commitment

Given \\(f\\) with coefficients \\(c_i\\), compute the commitment:

\\[
    C = \Pi_{i=0}^{d} (g_1^{(\alpha^{i})})^{c_i} = g_1^{\sum_{i=0}^{d} c_i \alpha^{i}} = g_1^{f(\alpha)}
\\]

Because the SRS contains \\(g_1^{(\alpha^{i})}\\), the prover can compute this multi-exponentiation without knowing \\(\alpha\\).

```rust
pub fn commit_kzg(p: &Polynomial<FqOrder>, pk: &PublicKeyKZG) -> CommitmentKZG {
    p.eval_with_powers_on_curve(&pk.alpha_1)
}
```

## Open

We want to open this polynomial at \\(u\\). Note that \\(f(x) - f(u)\\) has \\(x - u\\) as its root. Then we denote 

\begin{align*}
    f(x) - f(u) = (x - u) f_u(x)
\end{align*}

Since \\(f(x)\\) is \\(d\\)-th degree polynomial, \\(f_u\\) is a \\(d-1\\)-th degree polynomial.

Let 

\begin{align*}
    f_u(x) = \sum_{i=0}^{d-1} c' x^{i}
\end{align*}

Then, the witness is \\(w = \Pi_{i=0}^{d-1} (g_1^{(\alpha^{i})})^{c_{i}'}\\), which corresponds to \\(g_1^{\sum_{i=0}^{d-1} c'_i \alpha^{i}} = g_1^{f_u(\alpha)}\\)

## Verification

We want to verify that \\(f(\alpha) - f(u) = (\alpha - u)f_u(\alpha)\\). To check it, 

\begin{align*}
    e(c, g_2) / e(g_1, g_2)^{f(u)} \overset{?}{=} e(w, g_2^{\alpha} \cdot g_2^{-u})
\end{align*}

\\(e(c, g_2)\\) is equal to \\(e(g_1^{f(\alpha)}, g_2) = g_T^{f(\alpha)}\\), \\(e(g_1, g_2)^{f(u)} = g_T^{f(u)}\\) and \\(e(w, g_2^{\alpha} \cdot g_2^{-u})\\) is equal to \\(e(g_1^{f_u(\alpha)}, g_2^{\alpha-u}) = g_T^{f_u(\alpha)\cdot (\alpha - u)}\\).