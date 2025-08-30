# KZG

[Kate, Aniket, Gregory M. Zaverucha, and Ian Goldberg. "Constant-size commitments to polynomials and their applications." International conference on the theory and application of cryptology and information security. Berlin, Heidelberg: Springer Berlin Heidelberg, 2010.](https://link.springer.com/chapter/10.1007/978-3-642-17373-8_11)

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

Given \\(f\\) with coefficients \\(\\{c_i\\}_{i=0}^{d}\\), compute the commitment:

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

This \\(f_u\\) is a polynomial of degree \\(\leq d - 1\\) (because \\(X - u\\) divides \\(f(X) - f(u)\\)). Let \\(f_u(x) = \sum_{i=0}^{d-1} c'_{i} x^{i}\\). The prover then forms the witness:

\begin{align*}
    W = \Pi_{i=0}^{d-1} (g_1^{(\alpha^{i})})^{c_{i}'} = g_1^{\sum_{i=0}^{d-1} c'_i \alpha^{i}} = g_1^{f_u(\alpha)}
\end{align*}

Finally, the prover submits the evaluation and the witness \\(\langle f(u), W \rangle \\) to the verifier.

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

Here, the prover can computes \\(W\\) only using the public key, and \\(\alpha\\) is not required.

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

## Why this works:

### Correctness

Expand both sides in the target group \\(G_T\\):

- \\(e(c, g_2) = e(g_1^{f(\alpha)}, g_2) = g_T^{f(\alpha)}\\)
- \\(e(g_1, g_2)^{y} = g_T^{y} = g_T^{f(u)}\\)
- \\(e(W, g_2^{\alpha} \cdot g_2^{-u}) = e(g_1^{f_u(\alpha)}, g_2^{\alpha-u}) = g_T^{f_u(\alpha)\cdot (\alpha - u)}\\)


So the equality is exactly the exponent identity:

\begin{align*}
    f(\alpha) - f(u) = (\alpha - u) f_u(\alpha)
\end{align*}

which holds by the definition of \\(f_u(\alpha)\\). Thus a correct witness passes verification.

### Binding

The binding property of the KZG commitment protocol follows from the **t-Strong Diffie-Hellman (t-SDH)** assumption:

---

*Let \\(\alpha\\) be a random point in \\(\mathbb{F}_p\\), and let \\(g\\) be a generator of a group \\(G\\). Given the \\(d + 1\\)-tuple \\(\langle G, g, g^{\alpha}, g^{(\alpha^2)}, \dots g^{(\alpha^t)}  \rangle \in G^{d+1}\\), no polynomial-time adversary \\(\mathcal{A}\\) can output a pair \\(\langle c, g^{\frac{1}{\alpha + c}} \rangle\\) for any value of \\(c \in \mathbb{F}_p \setminus \\{-\alpha\\}\\) with non-negligible probability.*

---

Assume, for contradiction, that there exists an adversary \\(\mathcal{A}\\) that breaks the binding property of the KZG commitment. That is, \\(\mathcal{A}\\) produces two distinct valid openings of the same commitment \\(C\\) at the same evaluation point \\(u\\): \\(\langle y, W \rangle\\) and \\(\langle y', W' \rangle\\) with \\(y \neq y'\\). 

We construct a reduction \\(\mathcal{B}\\) that uses \\(\mathcal{A}\\) to solve the t-SDH problem.

Specifically, \\(\mathcal{B}\\) is given the t-SDH instance \\(\langle G_1, g_1, g_1^{\alpha}, \dots, g_1^{(\alpha^{t})} \rangle\\), and it sets the KZG public key to \\(\mathrm{PK} = (g_1, g_1^{\alpha}, \dots, g_1^{(\alpha^{t})}, g_2, g^{\alpha}_2)\\) and gives this \\(\mathrm{PK}\\) to \\(\mathcal{A}\\). Then, \\(\mathcal{A}\\) outputs a commitment \\(C\\) and two valid openings \\(\langle y, W \rangle\\) and \\(\langle y', W' \rangle\\) for the same \\(u \in \mathbb{F}_p\\) such that \\(y \neq y'\\). Both openings satisfy:

\begin{align*}
    e(C, g_2) = e(W, g_2^{\alpha - u}) \cdot e(g_1, g_2)^{y} = e(W', g_2^{\alpha - u}) \cdot e(g_1, g_2)^{y'}
\end{align*}

Since \\(g_1\\) is a generator of \\(G_1\\), there exist \\(\psi\\) and \\(\psi'\\) such that \\(W = g_1^{\psi}\\) and \\(W' = g_1^{\psi'}\\). Substitute into the pairing equations:

\begin{align*}
                &e(W, g_2^{\alpha - u}) \cdot e(g_1, g_2)^{y} = e(W', g^{\alpha - u}) \cdot e(g_1, g_2)^{y'} \\\\
    \Rightarrow &e(g_1^{\psi}, g_2^{\alpha - u}) \cdot e(g_1, g_2)^{y} = e(g_1^{\psi'}, g^{\alpha - u}) \cdot e(g_1, g_2)^{y'} \\\\
    \Rightarrow &g_T^{\psi \cdot (\alpha - u)} \cdot g_T^{y} = g_T^{\psi' \cdot (\alpha - u)} \cdot g_T^{y'} \\\\
    \Rightarrow &\psi(\alpha - u) + y = \psi'(\alpha - u) + y' \\\\
    \Rightarrow &(\psi - \psi')(\alpha - u) = y' - y \\\\
    \Rightarrow &\frac{\psi - \psi'}{y' - y} = \frac{1}{\alpha - u}
\end{align*}

Therefore, \\(\mathcal{B}\\) can compute

\begin{align*}
    (\frac{W}{W'})^{\frac{1}{y' - y}} = (g_1^{\psi - \psi'})^{\frac{1}{y' - y}} = g_1^{\frac{1}{\alpha - u}}
\end{align*}

and outputs \\(\langle -u, g_1^{\frac{1}{\alpha - u}} \rangle\\), which is a valid solution to the t-SDH problem. 

This constructino shows that \\(\mathcal{B}\\) can use \\(\mathcal{A}\\) to solve the t-SDH problem with the same success probability and withonly a constant-factor overhead in running time. Therefore, under the t-SDH assumption, the KZG commitment scheme is binding.

### Hiding

TBD

## Batch Proof

When the prover wants to open the same polynomial at multiple points, we can efficiently generate a witness that can prove those batched openings. Suppose the prover wants to open at \\(k\\) points \\((u_1, u_2, \dots, u_k)\\) whose evaluations are \\((y_1, y_2, \dots, y_k)\\). The setup and the commitment are the same with the above normal KZG commitment.

We then introduce two new polynomials for the opening process. Let \\(I(X)\\) be a polynomial that interpolates \\([(u_1, y_1), (u_2, y_2), \dots, (u_k, y_k)]\\), and \\(Z(X)\\) be \\(\prod_{i=1}^{k} (X - u_i)\\). Then, we construct the quotient polynomial as follows:

\begin{align*}
    f_u(X) = \frac{f(X) - I(X)}{Z(X)}
\end{align*}

Like the normal KZG commitment for a single opning point, the witness of the batch proof is defined as

\begin{align*}
    W = g_1^{f_u(\alpha)}
\end{align*}

Finally, the verifier checks the correctness as follows:

\begin{align*}
    e(W, g_2^{Z(\alpha)}) \overset{?}{=} e(C / I(\alpha), g_2)
\end{align*}

Note that we need to additionally include \\(\\{g_2^{(\alpha^i)}\\}_{i=2}^{k}\\) in the public key during the setup phase.

## Degree Bound Proof

When we want to ensure that the polunomial \\(f(X)\\) has degree at most \\(d\\), we can verify this property with a small modification to the standard KZG commitment scheme.

The idea is to define

\begin{align*}
    h(X) = X^{D - d} \cdot f(X)
\end{align*}

, where \\(D\\) is the maximum supported degree in the public key. The degree bound proof works as follows:

1. Along with the stanrad commitment \\(C\\), the prover also submits \\(g_1^{h(\alpha)}\\) to the verifier.
2. The verifier checks the pairing equation:

\begin{align*}
e(g_1^{h(\alpha)}, g_2) \overset{?}{=} e(C, g_2^{\alpha^{D - d}})
\end{align*}

If the prover is honest, then the equality holds because \\(e(g_1^{h(\alpha)}, g_2) = e(g_1, g_2)^{\alpha^{D - d} \cdot f(\alpha)}\\) and \\(e(C, g_2^{\alpha^{D - d}}) = e(g_1, g_2)^{f(\alpha) \cdot \alpha^{D - d}}\\).

If the degree of \\(f(X)\\) exceeds \\(d\\), then \\(h(X)\\) has degree larger than \\(D\\). In that case, the prover **cannot** compute \\(g_1^{h(\alpha)}\\) since the public key only contains up to \\(D\\)-th powers of \\(g_1^{\alpha}\\).

Note that the public key must also include \\(\\{g_2^{(\alpha^i)}\\}_{i=2}^{D}\\).

```rust
pub fn prove_degree_bound(
    p: &Polynomial<FqOrder>,
    pk: &PublicKeyKZG,
    d: usize,
) -> ProofDegreeBound {
    let max_d = pk.powers_1.len() - 1;
    let mut q_coef = (0..(max_d + 1 - d))
        .map(|_| FqOrder::zero())
        .collect::<Vec<_>>();
    q_coef[max_d - d] = FqOrder::one();
    let q = Polynomial { coef: q_coef };
    let r = p * &q;
    r.eval_with_powers_on_curve(&pk.powers_1)
}

pub fn verify_degree_bound(
    c: &CommitmentKZG,
    proof: &ProofDegreeBound,
    pk: &PublicKeyKZG,
    d: usize,
) -> bool {
    let max_d = pk.powers_1.len() - 1;
    optimal_ate_pairing(proof, &pk.powers_2[0]) == optimal_ate_pairing(c, &pk.powers_2[max_d - d])
}
```

See the following papers for the detail about the degree bound proof:

- Kohrita, Tohru, and Patrick Towa. "Zeromorph: Zero-knowledge multilinear-evaluation proofs from
homomorphic univariate commitments." Cryptology ePrint Archive (2023). https://eprint.iacr.org/2023/917
- Chiesa, Alessandro, Yuncong Hu, Mary Maller, et al. "Marlin: Preprocessing zkSNARKs with Universal
and Updatable SRS." Cryptology ePrint Archive (2019). https://eprint.iacr.org/2019/1047