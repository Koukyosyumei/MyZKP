# KZG

Let \\(G_1\\), \\(G_2\\), and \\(G_T\\) be groups such that there exists a bilinear mapping \\(e: G_1 \times G_2 \to G_T\\). Let \\(g_1 \in G_1\\) and \\(g_2 \in G_2\\) be fixed generators.

We want a short commitment to a polynomial over \\(\mathbb{F}\\) denoted as \\(f(X) = \sum_{i=0}^{d} c_{i} x^{i} \in \mathbb{F}_{q}[X]\\) with a public maximum degree bound \\(D\\) such that \\(d \leq D\\). The prover produces a commitment \\(C\\) for this polynomial so that it cannot change the polynomial after that. Then, the prover submits a short witness that proves the value \\(f(u)\\) at a point \\(u \in F_q\\)

## Setup

A trusted party samples a secret \\(\alpha \xleftarrow{\text{random}} \mathbb{F}_{q}\\), and publishes a public key \\(\mathrm{PK}\\):

\begin{align*}
    &\mathrm{PK} := (\\{g_1^{(\alpha^{i})}\\}^{D}_{i=0}, g_2, g^{\alpha}_2)
\end{align*}

The secret key \\(\alpha\\) should be securely destroyed after the generation of this \\(\mathrm{PK}\\).

```rust
pub struct PublicKeyKZG {
    pub powers_1: Vec<G1Point>,
    pub powers_2: Vec<G2Point>,
}

pub fn setup_kzg(g1: &G1Point, g2: &G2Point, n: usize) -> PublicKeyKZG {
    let alpha = FqOrder::random_element(&[]);

    let mut powers_1 = Vec::with_capacity(n);
    let mut alpha_power = FqOrder::one();
    for _ in 0..1 + n {
        powers_1.push(g1.mul_ref(alpha_power.clone().get_value()));
        alpha_power = alpha_power * alpha.clone();
    }

    let powers_2 = vec![g2.clone(), g2.mul_ref(alpha.get_value())];

    PublicKeyKZG { powers_1, powers_2 }
}
```

## Commitment

Given \\(f\\) with coefficients \\(c_i\\), compute the commitment:

\\[
    C = \Pi_{i=0}^{d} (g_1^{(\alpha^{i})})^{c_i} = g_1^{\sum_{i=0}^{d} c_i \alpha^{i}} = g_1^{f(\alpha)}
\\]

Because the public key contains \\(g_1^{(\alpha^{i})}\\), the prover can compute this multi-exponentiation without knowing \\(\alpha\\).

```rust
pub type CommitmentKZG = G1Point;

pub fn commit_kzg(p: &Polynomial<FqOrder>, pk: &PublicKeyKZG) -> CommitmentKZG {
    p.eval_with_powers_on_curve(&pk.powers_1)
}
```

## Open

To open \\(f\\) at \\(u\\), define the quotient polynomial:

\begin{align*}
    f_u(X) = \frac{f(X) - f(u)}{X - u}
\end{align*}

This \\(f_u\\) is a polynomial of degree \\(\leq d - 1\\) (because \\(X - u\\) divides \\(f(X) - f(u)\\)). Let \\(f_u(x) = \sum_{i=0}^{d-1} c' x^{i}\\). The prover then forms the witness:

\begin{align*}
    W = \Pi_{i=0}^{d-1} (g_1^{(\alpha^{i})})^{c_{i}'} = g_1^{\sum_{i=0}^{d-1} c'_i \alpha^{i}} = g_1^{f_u(\alpha)}
\end{align*}

Here, the prover can computes \\(W\\) only using the public key, and \\(\alpha\\) is not required.

```rust
pub struct ProofKZG {
    pub y: FqOrder,
    pub w: G1Point,
}

pub fn open_kzg(p: &Polynomial<FqOrder>, z: &FqOrder, pk: &PublicKeyKZG) -> ProofKZG {
    let y = p.eval(z);
    let y_poly = Polynomial {
        coef: (&[y.clone()]).to_vec(),
    };

    let q = (p - &y_poly) / Polynomial::<FqOrder>::from_monomials(&[z.clone()]);
    ProofKZG {
        y: y,
        w: q.eval_with_powers_on_curve(&pk.powers_1),
    }
}
```

## Verification

The verifier is given the public key, the commitment \\(C\\), the claimed evauation \\(y = f(u)\\), and the witness \\(W\\). They check the pairing equation:

\begin{align*}
    \frac{e(c, g_2)}{e(g_1, g_2)^{y}} \overset{?}{=} e(W, g_2^{\alpha} \cdot g_2^{-u})
\end{align*}

or equivalently:

\begin{align*}
    e(c, g_2) \overset{?}{=} e(g_1, g_2)^{y} \cdot e(W, g_2^{\alpha} \cdot g_2^{-u})
\end{align*}

```rust
pub fn verify_kzg(z: &FqOrder, c: &CommitmentKZG, proof: &ProofKZG, pk: &PublicKeyKZG) -> bool {
    let g1 = &pk.powers_1[0];
    let g2 = &pk.powers_2[0];
    let g2_s = &pk.powers_2[1];
    let g2_z = g2.mul_ref(z.clone().get_value());
    let g2_s_minus_z = g2_s.clone() - g2_z;

    let e1 = optimal_ate_pairing(&proof.w, &g2_s_minus_z);
    let e2 = optimal_ate_pairing(&g1, &g2);
    let e3 = optimal_ate_pairing(&c, &g2);

    e3 == e1 * (e2.pow(proof.y.get_value()))
}
```

### Why this works:

#### Correctness

Expand both sides in the target group \\(G_T\\):

- \\(e(c, g_2) = e(g_1^{f(\alpha)}, g_2) = g_T^{f(\alpha)}\\)
- \\(e(g_1, g_2)^{y} = g_T^{y} = g_T^{f(u)}\\)
- \\(e(W, g_2^{\alpha} \cdot g_2^{-u}) = e(g_1^{f_u(\alpha)}, g_2^{\alpha-u}) = g_T^{f_u(\alpha)\cdot (\alpha - u)}\\)


So the equality is exactly the exponent identity:

\begin{align*}
    f(\alpha) - f(u) = (\alpha - u) f_u(\alpha)
\end{align*}

which holds by the definition of \\(f_u(\alpha)\\). Thus a correct witness passes verification.

#### Binding

TBD

#### Hiding

TBD