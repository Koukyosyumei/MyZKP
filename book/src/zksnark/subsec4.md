# Bringing It All Together: SNARK

Let's recap the previous sections. First, the relationship between the inputs and outputs of any program can be expressed as a rank-one constraint system (R1CS) as follows:

\\[
  (L \cdot v) \circ (R \cdot v) - O \cdot v = 0  
\\]

, where \\(v\\) is the concatenation of all inputs, outputs, and intermediate values. This allows us to transform the statement, "I know the input values \\(x\\) that make the program returns the output values \\(y\\)", into "I know \\(v\\), whose outputs components are \\(y\\), that satisfies the constraint system corresponding to the program". 

Then, instead of separately checking each constraint (which corresponds to a row in the R1CS matrix), we can convert this into a more efficient polynomial-equivalence test.

## First Protocol: Naive Approach

The simplest protocol, based on the previous chapter, is as follows:

**Protocol (Setup)**

- **Interpolated Polynomial:** Construct \\(\\{\ell_i, r_i, o_i\\}_{i\in[d]}\\) from \\(L\\), \\(R\\), and \\(O\\), respectively.
- **Target Polynomial:** \\(t(x) = (x-1)(x-2) \cdots (x-m)\\)
- **Secret Seed:** A trusted party generates the random value \\(s\\) and \\(\alpha\\).
- **Proof Key:** Provided to the prover
  - \\(\\{g^{\ell_i(s)},g^{r_i(s)},g^{o_i(s)}\\}_{i\in[d]}\\)
  - \\(\\{g^{\alpha \ell_i(s)},g^{\alpha r_i(s)},g^{\alpha o_i(s)}\\}_{i\in[d]}\\)
  - \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Verification Key:**
  - \\(g^{t(s)}, g^{\alpha}\\)
- After distribution, the original secret seeds are securely destroyed.

Both the proof key and the verification key are publicly available, enabling anyone to generate and verify proofs based on the target program.

**Protocol (Proving)**

- Run the program to obtain the assignment vector \\(v\\).
- Compute the linear-combinations of polynomials
  - \\(\ell(x) = \sum_{i=1}^{d} v_i \ell_{i}(x),\quad r(x) = \sum_{i=1}^{d} v_i r_{i}(x),\quad o(x) = \sum_{i=1}^{d} v_i o_{i}(x)\\)
- Compute the quotient polynomial:
  - \\(h(x) = \frac{\ell(x) r(x) - o(x)}{t(x)}\\)
- Evaluate each polynomial at \\(s\\).
  - \\(g^{\ell(s)} = \prod^{d}_{i=1} (g^{\ell_i(s)})^{v _i} ,\quad g^{r(s)} = \prod^{d} _{i=1} (g^{r_i(s)})^{v_i} ,\quad g^{o(s)} = \prod^{d} _{i=1} (g^{o _i(s)})^{v _i} \\)
- Evaluate the shifted polynomials at \\(s\\).
  - \\(g^{\alpha \ell(s)} = \prod^{d} _{i=1} (g^{\alpha \ell _i(s)})^{v _i} ,\quad g^{\alpha r(s)} = \prod^{d} _{i=1} (g^{\alpha r _i(s)})^{v _i} ,\quad g^{\alpha o(s)} = \prod^{d} _{i=1} (g^{\alpha o_i(s)})^{v _i} \\)
- Compute \\(g^{h(s)}\\) from \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Proof**: 
  - \\((g^{\ell(s)}, g^{r(s)}, g^{o(s)}, g^{\alpha \ell(s)}, g^{\alpha r(s)}, g^{\alpha o(s)}, g^{h(s)})\\)

**Protocol (Verification)**

- Parse the proof as \\((g^{\ell}, g^r, g^o, g^{\ell'}, g^{r'}, g^{o'}, g^{h})\\)
- Check the polynomial restrictions
  - \\(e(g^{\ell}, g^{\alpha}) = e(g^{\ell'}, g),\quad e(g^{r}, g^{\alpha}) = e(g^{r'}, g),\quad e(g^{o}, g^{\alpha}) = e(g^{o'}, g)\\)
- Verify validity of the proof
  - \\(e(g^{\ell}, g^{r}) = e(g^t,g^h) \cdot e(g^o, g)\\)

**Implementation**

```rust
#[derive(Debug, Clone)]
pub struct ProofKey1 {
    g1_ell_i_vec: Vec<G1Point>,
    g1_r_i_vec: Vec<G1Point>,
    g2_r_i_vec: Vec<G2Point>,
    g1_o_i_vec: Vec<G1Point>,
    g1_alpha_ell_i_vec: Vec<G1Point>,
    g1_alpha_r_i_vec: Vec<G1Point>,
    g1_alpha_o_i_vec: Vec<G1Point>,
    g1_sj_vec: Vec<G1Point>,
}

#[derive(Debug, Clone)]
pub struct VerificationKey1 {
    g2_alpha: G2Point,
    g2_t_s: G2Point,
}

#[derive(Debug, Clone)]
pub struct Proof1 {
    g1_ell: G1Point,
    g1_r: G1Point,
    g2_r: G2Point,
    g1_o: G1Point,
    g1_ell_prime: G1Point,
    g1_r_prime: G1Point,
    g1_o_prime: G1Point,
    g1_h: G1Point,
}

pub fn setup(g1: &G1Point, g2: &G2Point, qap: &QAP<FqOrder>) -> (ProofKey1, VerificationKey1) {
    let s = FqOrder::random_element(&[]);
    let alpha = FqOrder::random_element(&[]);

    (
        ProofKey1 {
            g1_ell_i_vec: generate_challenge_vec(g1, &qap.ell_i_vec, &s),
            g1_r_i_vec: generate_challenge_vec(g1, &qap.r_i_vec, &s),
            g2_r_i_vec: generate_challenge_vec(g2, &qap.r_i_vec, &s),
            g1_o_i_vec: generate_challenge_vec(g1, &qap.o_i_vec, &s),
            g1_alpha_ell_i_vec: generate_alpha_challenge_vec(g1, &qap.ell_i_vec, &s, &alpha),
            g1_alpha_r_i_vec: generate_alpha_challenge_vec(g1, &qap.r_i_vec, &s, &alpha),
            g1_alpha_o_i_vec: generate_alpha_challenge_vec(g1, &qap.o_i_vec, &s, &alpha),
            g1_sj_vec: generate_s_powers(g1, &s, qap.m),
        },
        VerificationKey1 {
            g2_alpha: g2.mul_ref(alpha.get_value()),
            g2_t_s: g2.mul_ref(qap.t.eval(&s).sanitize().get_value()),
        },
    )
}

pub fn prove(assignment: &Vec<FqOrder>, proof_key: &ProofKey1, qap: &QAP<FqOrder>) -> Proof1 {
    Proof1 {
        g1_ell: accumulate_curve_points(&proof_key.g1_ell_i_vec, assignment),
        g1_r: accumulate_curve_points(&proof_key.g1_r_i_vec, assignment),
        g2_r: accumulate_curve_points(&proof_key.g2_r_i_vec, assignment),
        g1_o: accumulate_curve_points(&proof_key.g1_o_i_vec, assignment),
        g1_ell_prime: accumulate_curve_points(&proof_key.g1_alpha_ell_i_vec, assignment),
        g1_r_prime: accumulate_curve_points(&proof_key.g1_alpha_r_i_vec, assignment),
        g1_o_prime: accumulate_curve_points(&proof_key.g1_alpha_o_i_vec, assignment),
        g1_h: get_h(qap, assignment).eval_with_powers_on_curve(&proof_key.g1_sj_vec),
    }
}

pub fn verify(
    g1: &G1Point,
    g2: &G2Point,
    proof: &Proof1,
    verification_key: &VerificationKey1,
) -> bool {
    let pairing1 = optimal_ate_pairing(&proof.g1_ell, &verification_key.g2_alpha);
    let pairing2 = optimal_ate_pairing(&proof.g1_ell_prime, &g2);
    if pairing1 != pairing2 {
        return false;
    }

    let pairing3 = optimal_ate_pairing(&proof.g1_r, &verification_key.g2_alpha);
    let pairing4 = optimal_ate_pairing(&proof.g1_r_prime, &g2);
    if pairing3 != pairing4 {
        return false;
    }

    let pairing5 = optimal_ate_pairing(&proof.g1_o, &verification_key.g2_alpha);
    let pairing6 = optimal_ate_pairing(&proof.g1_o_prime, &g2);
    if pairing5 != pairing6 {
        return false;
    }

    let pairing7 = optimal_ate_pairing(&proof.g1_ell, &proof.g2_r);
    let pairing8 = optimal_ate_pairing(&proof.g1_h, &verification_key.g2_t_s);
    let pairing9 = optimal_ate_pairing(&proof.g1_o, &g2);

    pairing7 == pairing8 * pairing9
}

pub fn generate_challenge_vec<F1: Field, F2: Field, E: EllipticCurve>(
    point: &EllipticCurvePoint<F1, E>,
    poly_vec: &[Polynomial<F2>],
    s: &F2,
) -> Vec<EllipticCurvePoint<F1, E>> {
    poly_vec
        .iter()
        .map(|poly| point.mul_ref(poly.eval(s).sanitize().get_value()))
        .collect()
}

pub fn generate_alpha_challenge_vec<F1: Field, F2: Field, E: EllipticCurve>(
    point: &EllipticCurvePoint<F1, E>,
    poly_vec: &[Polynomial<F2>],
    s: &F2,
    alpha: &F2,
) -> Vec<EllipticCurvePoint<F1, E>> {
    poly_vec
        .iter()
        .map(|poly| point.mul_ref((alpha.mul_ref(&poly.eval(s).sanitize())).get_value()))
        .collect()
}

pub fn generate_s_powers<F1: Field, F2: Field, E: EllipticCurve>(
    point: &EllipticCurvePoint<F1, E>,
    s: &F2,
    m: usize,
) -> Vec<EllipticCurvePoint<F1, E>> {
    let mut powers = Vec::with_capacity(m + 1);
    let mut current = F2::one();
    for _ in 0..=m {
        powers.push(point.mul_ref(current.get_value()));
        current = current * s.clone();
    }
    powers
}

pub fn accumulate_curve_points<F1: Field, F2: Field, E: EllipticCurve>(
    g_vec: &[EllipticCurvePoint<F1, E>],
    assignment: &[F2],
) -> EllipticCurvePoint<F1, E> {
    g_vec.iter().zip(assignment.iter()).fold(
        EllipticCurvePoint::<F1, E>::point_at_infinity(),
        |acc, (g, &ref a)| acc + g.mul_ref(a.get_value()),
    )
}

pub fn accumulate_polynomials<F: Field>(
    poly_vec: &[Polynomial<F>],
    assignment: &[F],
) -> Polynomial<F> {
    poly_vec
        .iter()
        .zip(assignment.iter())
        .fold(Polynomial::<F>::zero(), |acc, (poly, &ref a)| {
            acc + poly.clone() * a.clone()
        })
}

pub fn get_h<F: Field>(qap: &QAP<F>, assignment: &Vec<F>) -> Polynomial<F> {
    let ell = accumulate_polynomials(&qap.ell_i_vec, assignment);
    let r = accumulate_polynomials(&qap.r_i_vec, assignment);
    let o = accumulate_polynomials(&qap.o_i_vec, assignment);
    (ell * r - o) / qap.t.clone()
}
```

**Vulnerability**

A critical issue with this protocol is that the checks may pass even if \\(g^{\ell}\\) is computed not from \\(\\{g^{\ell_i(s)}\\}_{i \in [d]}\\) but from \\(\\{g^{r_i(s)}\\} _{i \in [d]}\\), \\(\\{g^{o_i(s)}\\} _{i \in [d]}\\), or their combinations. The same issue applies to \\(g^{r}\\) and \\(g^{o}\\). 

For example, if the prover sends \\((g^{\ell(s)}, g^{\ell(s)}, g^{o(s)}, g^{\alpha r(s)}, g^{\alpha \ell(s)}, g^{\alpha o(s)}, g^{h(s)})\\) as the proof, all the verification checks still pass, even although the proved statement differs from the original one.

```rust
pub fn interchange_attack(proof: &Proof1) -> Proof1 {
    let mut new_proof = proof.clone();
    new_proof.g1_r = proof.g1_ell.clone();
    new_proof.g1_r_prime = proof.g1_ell_prime.clone();
    new_proof
}
```

## Second Protocol: Non-Interchangibility

To address the interchangeability issue, the next protocol uses distinct the different \\(\alpha\\)-shift for \\(\ell\\), \\(r\\), and \\(o\\).

**Protocol (Setup)**

- **Interpolated Polynomial:** Construct \\(\\{\ell_i, r_i, o_i\\}_{i\in[d]}\\) from \\(L\\), \\(R\\), and \\(O\\), respectively.
- **Target Polynomial:** \\(t(x) = (x-1)(x-2) \cdots (x-m)\\)
- **Secret Seed:** A trusted party generates the random value \\(s\\), *\\(\alpha_{\ell}\\), \\(\alpha_r\\), and \\(\alpha_o\\)*.
- **Proof Key (for the prover):**
  - \\(\\{g^{\ell_i(s)},g^{r_i(s)},g^{o_i(s)}\\}_{i\in[d]}\\)
  - *\\(\\{g^{\alpha_{\ell} \ell_i(s)},g^{\alpha_{r} r_i(s)},g^{\alpha_{o} o_i(s)}\\}_{i\in[d]}\\)*
  - \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Verification Key (public):**
  - \\(g^{t(s)}\\) *\\(, g^{\alpha_{\ell}},g^{\alpha_{r}},g^{\alpha_{o}}\\)*
- After distribution, the original secret seeds are securely destroyed.

**Protocol (Proving)**

- Execute the program and get the assignment \\(v\\).
- Compute the linear-combinations of polynomials
  - \\(\ell(x) = \sum_{i=1}^{d} v_i \ell_{i}(x),\quad r(x) = \sum_{i=1}^{d} v_i r_{i}(x),\quad o(x) = \sum_{i=1}^{d} v_i o_{i}(x)\\)
- Compute the quotient polynomial:
  - \\(h(x) = \frac{\ell(x) r(x) - o(x)}{t(x)}\\)
- Evaluate each polynomial at \\(s\\).
  - \\(g^{\ell(s)} = \prod^{d}_{i=1} (g^{\ell_i(s)})^{v _i} ,\quad g^{r(s)} = \prod^{d} _{i=1} (g^{r_i(s)})^{v_i} ,\quad g^{o(s)} = \prod^{d} _{i=1} (g^{o _i(s)})^{v _i} \\)
- Evaluate each shifted polynomial at \\(s\\).
  - *\\(g^{\alpha_{\ell} \ell(s)} = \prod^{d}_{i=1} (g^{\alpha _{\ell} \ell_i(s)})^{v_i} \\)*, *\\(g^{\alpha_{r} r(s)} = \prod^{d}_{i=1} (g^{\alpha _{r} r_i(s)})^{v_i} \\)*, *\\(g^{\alpha_{o} o(s)} = \prod^{d}_{i=1} (g^{\alpha _{o} o_i(s)})^{v_i} \\)*
- Calculate \\(g^{h(s)}\\) from \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Proof**: 
  - \\((g^{\ell(s)}, g^{r(s)}, g^{o(s)},\\) *\\(g^{\alpha_{\ell} \ell(s)}, g^{\alpha_{r} r(s)}, g^{\alpha_{o} o(s)},\\)* \\(g^{h(s)})\\)

**Protocol (Verification)**

- Parse proof as \\((g^{\ell}, g^r, g^o, g^{\ell'}, g^{r'}, g^{o'}, g^{h})\\)
- Check polynomial restrictions
  - *\\(e(g^{\ell}, g^{\alpha_{\ell}}) = e(g^{\ell'}, g)\\)*, \\(\quad \\) *\\(e(g^{r}, g^{\alpha_{r}}) = e(g^{r'}, g)\\)*, \\(\quad \\) *\\(e(g^{o}, g^{\alpha_{o}}) = e(g^{o'}, g)\\)*
- Validity check
  - \\(e(g^{\ell}, g^{r}) = e(g^t,g^h) \cdot e(g^o, g)\\)

**Implementation**

```rust
#[derive(Debug, Clone)]
pub struct ProofKey2 {
    g1_ell_i_vec: Vec<G1Point>,
    g1_r_i_vec: Vec<G1Point>,
    g2_r_i_vec: Vec<G2Point>,
    g1_o_i_vec: Vec<G1Point>,
    g1_alpha_ell_i_vec: Vec<G1Point>,
    g1_alpha_r_i_vec: Vec<G1Point>,
    g1_alpha_o_i_vec: Vec<G1Point>,
    g1_sj_vec: Vec<G1Point>,
}

#[derive(Debug, Clone)]
pub struct VerificationKey2 {
    g2_alpha_ell: G2Point,
    g2_alpha_r: G2Point,
    g2_alpha_o: G2Point,
    g2_t_s: G2Point,
}

#[derive(Debug, Clone)]
pub struct Proof2 {
    g1_ell: G1Point,
    g1_r: G1Point,
    g2_r: G2Point,
    g1_o: G1Point,
    g1_ell_prime: G1Point,
    g1_r_prime: G1Point,
    g1_o_prime: G1Point,
    g1_h: G1Point,
}

pub fn setup(g1: &G1Point, g2: &G2Point, qap: &QAP<FqOrder>) -> (ProofKey2, VerificationKey2) {
    let s = FqOrder::random_element(&[]);
    let alpha_ell = FqOrder::random_element(&[]);
    let alpha_r = FqOrder::random_element(&[]);
    let alpha_o = FqOrder::random_element(&[]);

    (
        ProofKey2 {
            g1_ell_i_vec: generate_challenge_vec(g1, &qap.ell_i_vec, &s),
            g1_r_i_vec: generate_challenge_vec(g1, &qap.r_i_vec, &s),
            g2_r_i_vec: generate_challenge_vec(g2, &qap.r_i_vec, &s),
            g1_o_i_vec: generate_challenge_vec(g1, &qap.o_i_vec, &s),
            g1_alpha_ell_i_vec: generate_alpha_challenge_vec(g1, &qap.ell_i_vec, &s, &alpha_ell),
            g1_alpha_r_i_vec: generate_alpha_challenge_vec(g1, &qap.r_i_vec, &s, &alpha_r),
            g1_alpha_o_i_vec: generate_alpha_challenge_vec(g1, &qap.o_i_vec, &s, &alpha_o),
            g1_sj_vec: generate_s_powers(g1, &s, qap.m),
        },
        VerificationKey2 {
            g2_alpha_ell: g2.mul_ref(alpha_ell.get_value()),
            g2_alpha_r: g2.mul_ref(alpha_r.get_value()),
            g2_alpha_o: g2.mul_ref(alpha_o.get_value()),
            g2_t_s: g2.mul_ref(qap.t.eval(&s).sanitize().get_value()),
        },
    )
}

pub fn prove(assignment: &Vec<FqOrder>, proof_key: &ProofKey2, qap: &QAP<FqOrder>) -> Proof2 {
    Proof2 {
        g1_ell: accumulate_curve_points(&proof_key.g1_ell_i_vec, assignment),
        g1_r: accumulate_curve_points(&proof_key.g1_r_i_vec, assignment),
        g2_r: accumulate_curve_points(&proof_key.g2_r_i_vec, assignment),
        g1_o: accumulate_curve_points(&proof_key.g1_o_i_vec, assignment),
        g1_ell_prime: accumulate_curve_points(&proof_key.g1_alpha_ell_i_vec, assignment),
        g1_r_prime: accumulate_curve_points(&proof_key.g1_alpha_r_i_vec, assignment),
        g1_o_prime: accumulate_curve_points(&proof_key.g1_alpha_o_i_vec, assignment),
        g1_h: get_h(qap, assignment).eval_with_powers_on_curve(&proof_key.g1_sj_vec),
    }
}

pub fn verify(
    g1: &G1Point,
    g2: &G2Point,
    proof: &Proof2,
    verification_key: &VerificationKey2,
) -> bool {
    let pairing1 = optimal_ate_pairing(&proof.g1_ell, &verification_key.g2_alpha_ell);
    let pairing2 = optimal_ate_pairing(&proof.g1_ell_prime, &g2);
    if pairing1 != pairing2 {
        return false;
    }

    let pairing3 = optimal_ate_pairing(&proof.g1_r, &verification_key.g2_alpha_r);
    let pairing4 = optimal_ate_pairing(&proof.g1_r_prime, &g2);
    if pairing3 != pairing4 {
        return false;
    }

    let pairing5 = optimal_ate_pairing(&proof.g1_o, &verification_key.g2_alpha_o);
    let pairing6 = optimal_ate_pairing(&proof.g1_o_prime, &g2);
    if pairing5 != pairing6 {
        return false;
    }

    let pairing7 = optimal_ate_pairing(&proof.g1_ell, &proof.g2_r);
    let pairing8 = optimal_ate_pairing(&proof.g1_h, &verification_key.g2_t_s);
    let pairing9 = optimal_ate_pairing(&proof.g1_o, &g2);

    pairing7 == pairing8 * pairing9
}
```

**Vulnerability**

This protocol resolves interchangeability but does not enforce consistency acros \\(\ell\\), \\(r\\), and \\(o\\). Variables \\(v_i\\) can still take different values in each polynoimal because verification checks are performed separately.

```rust
pub fn inconsistent_variable_attack(
    assignment_ell: &Vec<FqOrder>,
    assignment_r: &Vec<FqOrder>,
    assignment_o: &Vec<FqOrder>,
    proof_key: &ProofKey2,
    qap: &QAP<FqOrder>,
) -> Proof2 {
    let ell = accumulate_polynomials(&qap.ell_i_vec, assignment_ell);
    let r = accumulate_polynomials(&qap.r_i_vec, assignment_r);
    let o = accumulate_polynomials(&qap.o_i_vec, assignment_o);
    let h = (ell * r - o) / qap.t.clone();

    Proof2 {
        g1_ell: accumulate_curve_points(&proof_key.g1_ell_i_vec, assignment_ell),
        g1_r: accumulate_curve_points(&proof_key.g1_r_i_vec, assignment_r),
        g2_r: accumulate_curve_points(&proof_key.g2_r_i_vec, assignment_r),
        g1_o: accumulate_curve_points(&proof_key.g1_o_i_vec, assignment_o),
        g1_ell_prime: accumulate_curve_points(&proof_key.g1_alpha_ell_i_vec, assignment_ell),
        g1_r_prime: accumulate_curve_points(&proof_key.g1_alpha_r_i_vec, assignment_r),
        g1_o_prime: accumulate_curve_points(&proof_key.g1_alpha_o_i_vec, assignment_o),
        g1_h: h.eval_with_powers_on_curve(&proof_key.g1_sj_vec),
    }
}
```

## Thrid Protocol: Variable Consistency

To achive the variable consistency, we employ a checksum mechanism. Specifically, we first draw a new random value \\(\beta\\) and define the checksum of \\(v_i\\) as \\(g^{\beta(\ell_{i}(s) + r_i(s) + o_i(s))}\\). Let \\(v _{\ell, i}\\), \\(v _{r,i}\\), \\(v _{o,i}\\), and \\(v _{\beta,i}\\) denote the \\(i\\)-th value of the assignment vectors for \\(\ell\\), \\(r\\), \\(o\\) and the checksum, respectively. If all of them are the same, the following equation holds:

\\[
e(g^{v _{\ell, i} \ell_i(s)} g^{v _{r, i} r_i(s)} g^{v _{o, i} o_i(s)}, g^{\beta}) = e(g^{v _{\beta, i} \beta(\ell _{i}(s) + r _{i}(s) + o _{i}(s))}, g)
\\]

Unfortunately, this condition is not strictly equivalent. For example, consider the case where \\(\ell_i(x) = r_i(x)\\). In this scenario, we have:

\begin{align*}
  &\beta(v _{\ell, i} \ell_i(s) + v _{r, i} r_i(s) + v _{o, i} o_i(s)) = v _{\beta, i} \beta (\ell _{i}(s) + r _{i}(s) + o _{i}(s)) \\\\
  \iff &\beta(v _{\ell, i} \ell_i(s) + v _{r, i} \ell_i(s) + v _{o, i} o_i(s)) = v _{\beta, i} \beta (2\ell _{i}(s) + o _{i}(s))
\end{align*}

This equation holds for arbitrary \\(v _{r,i}\\) and \\(v _{o,i}\\) if we set \\(v _{\beta, i} = v _{o, i}\\) and \\(v _{\ell, i} = 2 v _{o, i} - v _{r, i}\\).

To address this issue, we use distinct different \\(\beta\\) values for \\(\ell\\), \\(r\\) and \\(o\\). The consistency check then verifies the following equation:

\\[
e(g^{v _{\ell, i} \ell _{i}(s)}, g^{\beta _{\ell}}) \cdot e(g^{v _{r, i} r _{i}(s)}, g^{\beta _{r}}) \cdot e(g^{v _{o, i} o _{i}(s)}, g^{\beta _{o}}) = e(g^{v _{\beta, i}(\beta _{\ell} \ell _{i}(s) + \beta _{r} r _{i}(s) + \beta _{o} o _{i}(s))}, g)   
\\]

The new protocol using the above variable-consistency check is as follows:

**Protocol (Setup)**

- **Interpolated Polynomial:** Construct \\(\\{\ell_i, r_i, o_i\\}_{i\in[d]}\\) from \\(L\\), \\(R\\), and \\(O\\), respectively.
- **Target Polynomial:** \\(t(x) = (x-1)(x-2) \cdots (x-m)\\)
- **Secret Seed:** A trusted party generates the random value \\(s\\), \\(\alpha_{\ell}\\), \\(\alpha_r\\), \\(\alpha_o\\), *\\(\beta_{\ell}\\), \\(\beta_{r}\\), and \\(\beta_{o}\\)*.
- **Proof Key:** Provided to the prover
  - \\(\\{g^{\ell_i(s)},g^{r_i(s)},g^{o_i(s)}\\}_{i\in[d]}\\)
  - \\(\\{g^{\alpha_{\ell} \ell_i(s)},g^{\alpha_{r} r_i(s)},g^{\alpha_{o} o_i(s)}\\}_{i\in[d]}\\)
  - \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
  - *\\(\\{g^{\beta _{\ell} \ell _{i}(s) + \beta _{r} r _{i}(s) + \beta _{o} o _{i}(s)}\\} _{i \in [d]}\\)*
- **Verification Key:**
  - \\(g^{t(s)}, g^{\alpha_{\ell}},g^{\alpha_{r}},g^{\alpha_{o}}\\) *\\(, g^{\beta_{\ell}}, g^{\beta_{r}}, g^{\beta_{o}}\\)*
- After distribution, the original secret seeds are securely destroyed.

**Protocol (Proving)**

- Execute the program and get the assignment vector \\(v\\).
- Compute the linear-combinations of polynomials
  - \\(\ell(x) = \sum_{i=1}^{d} v_i \ell_{i}(x),\quad r(x) = \sum_{i=1}^{d} v_i r_{i}(x),\quad o(x) = \sum_{i=1}^{d} v_i o_{i}(x)\\)
- Compute the quotient polynomial:
  - \\(h(x) = \frac{\ell(x) r(x) - o(x)}{t(x)}\\)
- Evaluate each polynomial at \\(s\\):
  - \\(g^{\ell(s)} = \prod^{d}_{i=1} (g^{\ell_i(s)})^{v _i} ,\quad g^{r(s)} = \prod^{d} _{i=1} (g^{r_i(s)})^{v_i} ,\quad g^{o(s)} = \prod^{d} _{i=1} (g^{o _i(s)})^{v _i} \\)
- Evaluate each shifted polynomial at \\(s\\):
  - \\(g^{\alpha _{\ell} \ell(s)} = \prod^{d} _{i=1} (g^{\alpha _{\ell} \ell _i(s)})^{v _i} ,g^{\alpha _{r} r(s)} = \prod^{d} _{i=1} (g^{\alpha _{r} r _i(s)})^{v _i} ,g^{\alpha _{o} o(s)} = \prod^{d} _{i=1} (g^{\alpha _{o} o _i(s)})^{v _i} \\)
- *Evaluate each consistency polynomial at \\(s\\):*
  - *\\(g^{z(s)} = \prod^{d}_{i=1} (g^{\beta _{\ell} \ell _{i}(s) + \beta _{r} r _{i}(s) + \beta _{o} o _{i}(s)})^{v _{i}}\\)*
- Calculate \\(g^{h(s)}\\) from \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Proof**: 
  - \\((\\) \\(g^{\ell(s)}, g^{r(s)}, g^{o(s)}, g^{\alpha_{\ell} \ell(s)}, g^{\alpha_{r} r(s)}, g^{\alpha_{o} o(s)}, g^{h(s)},\\) *\\(g^{z(s)}\\)* \\()\\)

**Protocol (Verification)**

- Parse proof as \\((g^{\ell}, g^r, g^o, g^{\ell'}, g^{r'}, g^{o'}, g^{h}, g^{z})\\)
- Check polynomial restrictions
  - \\(\quad e(g^{\ell}, g^{\alpha_{\ell}}) = e(g^{\ell'}, g),\quad e(g^{r}, g^{\alpha_{r}}) = e(g^{r'}, g),\quad e(g^{o}, g^{\alpha_{o}}) = e(g^{o'}, g)\\)
- Check variable consistency
  - *\\(e(g^{\ell}, g^{\beta _{\ell}}) \cdot e(g^{r}, g^{\beta _{r}}) \cdot e(g^{o}, g^{\beta _{o}}) = e(g^{z}, g)\\)*
- Validity check
  - \\(e(g^{\ell}, g^{r}) = e(g^t,g^h) \cdot e(g^o, g)\\)

**Implementation**

```rust
#[derive(Debug, Clone)]
pub struct ProofKey3 {
    g1_ell_i_vec: Vec<G1Point>,
    g1_r_i_vec: Vec<G1Point>,
    g2_r_i_vec: Vec<G2Point>,
    g1_o_i_vec: Vec<G1Point>,
    g1_alpha_ell_i_vec: Vec<G1Point>,
    g1_alpha_r_i_vec: Vec<G1Point>,
    g1_alpha_o_i_vec: Vec<G1Point>,
    g1_sj_vec: Vec<G1Point>,
    g1_checksum_vec: Vec<G1Point>,
}

#[derive(Debug, Clone)]
pub struct VerificationKey3 {
    g2_alpha_ell: G2Point,
    g2_alpha_r: G2Point,
    g2_alpha_o: G2Point,
    g2_beta_ell: G2Point,
    g2_beta_r: G2Point,
    g2_beta_o: G2Point,
    g2_t_s: G2Point,
}

#[derive(Debug, Clone)]
pub struct Proof3 {
    g1_ell: G1Point,
    g1_r: G1Point,
    g2_r: G2Point,
    g1_o: G1Point,
    g1_ell_prime: G1Point,
    g1_r_prime: G1Point,
    g1_o_prime: G1Point,
    g1_h: G1Point,
    g1_z: G1Point,
}

pub fn setup(g1: &G1Point, g2: &G2Point, qap: &QAP<FqOrder>) -> (ProofKey3, VerificationKey3) {
    let s = FqOrder::random_element(&[]);
    let alpha_ell = FqOrder::random_element(&[]);
    let alpha_r = FqOrder::random_element(&[]);
    let alpha_o = FqOrder::random_element(&[]);
    let beta_ell = FqOrder::random_element(&[]);
    let beta_r = FqOrder::random_element(&[]);
    let beta_o = FqOrder::random_element(&[]);

    let mut g1_checksum_vec = Vec::with_capacity(qap.d);

    for i in 0..qap.d {
        let ell_i_s = qap.ell_i_vec[i].eval(&s).sanitize();
        let r_i_s = qap.r_i_vec[i].eval(&s).sanitize();
        let o_i_s = qap.o_i_vec[i].eval(&s).sanitize();
        g1_checksum_vec.push(
            g1.mul_ref(
                ((beta_ell.mul_ref(&ell_i_s))
                    .add_ref(&beta_r.mul_ref(&r_i_s))
                    .add_ref(&beta_o.mul_ref(&o_i_s)))
                .get_value(),
            ),
        );
    }

    (
        ProofKey3 {
            g1_ell_i_vec: generate_challenge_vec(g1, &qap.ell_i_vec, &s),
            g1_r_i_vec: generate_challenge_vec(g1, &qap.r_i_vec, &s),
            g2_r_i_vec: generate_challenge_vec(g2, &qap.r_i_vec, &s),
            g1_o_i_vec: generate_challenge_vec(g1, &qap.o_i_vec, &s),
            g1_alpha_ell_i_vec: generate_alpha_challenge_vec(g1, &qap.ell_i_vec, &s, &alpha_ell),
            g1_alpha_r_i_vec: generate_alpha_challenge_vec(g1, &qap.r_i_vec, &s, &alpha_r),
            g1_alpha_o_i_vec: generate_alpha_challenge_vec(g1, &qap.o_i_vec, &s, &alpha_o),
            g1_sj_vec: generate_s_powers(g1, &s, qap.m),
            g1_checksum_vec: g1_checksum_vec,
        },
        VerificationKey3 {
            g2_alpha_ell: g2.mul_ref(alpha_ell.get_value()),
            g2_alpha_r: g2.mul_ref(alpha_r.get_value()),
            g2_alpha_o: g2.mul_ref(alpha_o.get_value()),
            g2_beta_ell: g2.mul_ref(beta_ell.get_value()),
            g2_beta_r: g2.mul_ref(beta_r.get_value()),
            g2_beta_o: g2.mul_ref(beta_o.get_value()),
            g2_t_s: g2.mul_ref(qap.t.eval(&s).sanitize().get_value()),
        },
    )
}

pub fn prove(assignment: &Vec<FqOrder>, proof_key: &ProofKey3, qap: &QAP<FqOrder>) -> Proof3 {
    Proof3 {
        g1_ell: accumulate_curve_points(&proof_key.g1_ell_i_vec, assignment),
        g1_r: accumulate_curve_points(&proof_key.g1_r_i_vec, assignment),
        g2_r: accumulate_curve_points(&proof_key.g2_r_i_vec, assignment),
        g1_o: accumulate_curve_points(&proof_key.g1_o_i_vec, assignment),
        g1_ell_prime: accumulate_curve_points(&proof_key.g1_alpha_ell_i_vec, assignment),
        g1_r_prime: accumulate_curve_points(&proof_key.g1_alpha_r_i_vec, assignment),
        g1_o_prime: accumulate_curve_points(&proof_key.g1_alpha_o_i_vec, assignment),
        g1_h: get_h(qap, assignment).eval_with_powers_on_curve(&proof_key.g1_sj_vec),
        g1_z: accumulate_curve_points(&proof_key.g1_checksum_vec, assignment),
    }
}

pub fn verify(
    g1: &G1Point,
    g2: &G2Point,
    proof: &Proof3,
    verification_key: &VerificationKey3,
) -> bool {
    let pairing1 = optimal_ate_pairing(&proof.g1_ell, &verification_key.g2_alpha_ell);
    let pairing2 = optimal_ate_pairing(&proof.g1_ell_prime, &g2);
    if pairing1 != pairing2 {
        return false;
    }

    let pairing3 = optimal_ate_pairing(&proof.g1_r, &verification_key.g2_alpha_r);
    let pairing4 = optimal_ate_pairing(&proof.g1_r_prime, &g2);
    if pairing3 != pairing4 {
        return false;
    }

    let pairing5 = optimal_ate_pairing(&proof.g1_o, &verification_key.g2_alpha_o);
    let pairing6 = optimal_ate_pairing(&proof.g1_o_prime, &g2);
    if pairing5 != pairing6 {
        return false;
    }

    let pairing7 = optimal_ate_pairing(&proof.g1_ell, &proof.g2_r);
    let pairing8 = optimal_ate_pairing(&proof.g1_h, &verification_key.g2_t_s);
    let pairing9 = optimal_ate_pairing(&proof.g1_o, &g2);

    if pairing7 != pairing8 * pairing9 {
        return false;
    }

    let pairing10 = optimal_ate_pairing(&proof.g1_ell, &verification_key.g2_beta_ell);
    let pairing11 = optimal_ate_pairing(&proof.g1_r, &verification_key.g2_beta_r);
    let pairing12 = optimal_ate_pairing(&proof.g1_o, &verification_key.g2_beta_o);
    let pairing13 = optimal_ate_pairing(&proof.g1_z, &g2);

    pairing10 * pairing11 * pairing12 == pairing13
}
```

**Vulnerability**

Despite these checks, the protocol is vulnerable to **malleability**. Specifically, a malicious prover can exploit the polynomial restriction check to introduce arbitrary constants, altering the proof witout detection.

Recall that the verifier validates whether the submitted \\(g^{\ell}\\) is actually calculated by \\(\\{g^{\ell _i(s)}\\} _{i \in [d]}\\) by checking \\(e(g^{\ell}, g^{\alpha _{\ell}}) = e(g^{\ell'}, g)\\). However, this process is not sound. Recall that the verification key is publicaly available, and the prover knows both of \\(g^{\alpha _{\ell}}\\) and \\(g^{\beta _{\ell}}\\). Here, suppose the prover submits \\(g^{\ell} g^{c}\\) and \\(g^{\ell'} (g^{\alpha _{\ell}})^{c}\\) insteads of \\(g^{\ell}\\) and \\(g^{\ell'}\\), where \\(c\\) is a constatn value. Then, the polynomial restriction check still passes:

\\[
e(g^{\ell} g^{c}, g^{\alpha _{\ell}}) = e(g^{\alpha _{\ell} \ell + \alpha _{\ell} c}, g) = e(g^{\ell'}g^{\alpha _{\ell}c}, g) = e(g^{\ell'} (g^{\alpha _{\ell}})^{c}, g)
\\]

In addition, if the prover submits \\(g^{z} (g^{\beta _{\ell}})^{c}\\) as the checksum, it also passes the polynomial checksum verification:

\\[
e(g^{\ell} g^{c}, g^{\beta _{\ell}}) \cdot e(g^{r}, g^{\beta _{r}}) \cdot e(g^{o}, g^{\beta _{o}}) = e(g^{z} (g^{\beta _{\ell}})^{c}, g)
\\]

This phenomenon also can occur for \\(r\\) and \\(o\\).

## Forth Protocol: Non-Malleability

One way to surrogate the above malleability is hiding \\(g^{\beta _{\ell}}\\), \\(g^{\beta _{r}}\\), and \\(g^{\beta _{o}}\\) by powering them with a new random value \\(\eta\\).


**Protocol (Setup)**

- **Interpolated Polynomial:** Construct \\(\\{\ell_i, r_i, o_i\\}_{i\in[d]}\\) from \\(L\\), \\(R\\), and \\(O\\), respectively.
- **Target Polynomial:** \\(t(x) = (x-1)(x-2) \cdots (x-m)\\)
- **Secret Seed:** A trusted party generates the random value \\(s\\), \\(\alpha_{\ell}\\), \\(\alpha_r\\), \\(\alpha_o\\), \\(\beta_{\ell}\\), \\(\beta_{r}\\), \\(\beta_{o}\\), and *\\(\eta\\)*.
- **Proof Key:** Provided to the prover
  - \\(\\{g^{\ell_i(s)},g^{r_i(s)},g^{o_i(s)}\\}_{i\in[d]}\\)
  - \\(\\{g^{\alpha_{\ell} \ell_i(s)},g^{\alpha_{r} r_i(s)},g^{\alpha_{o} o_i(s)}\\}_{i\in[d]}\\)
  - \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
  - \\(\\{g^{\beta _{\ell} \ell _{i}(s) + \beta _{r} r _{i}(s) + \beta _{o} o _{i}(s)}\\} _{i \in [d]}\\)
- **Verification Key:**
  - \\(g^{t(s)}, g^{\alpha_{\ell}},g^{\alpha_{r}},g^{\alpha_{o}}\\) *\\(, g^{\beta_{\ell} \eta}, g^{\beta_{r} \eta}, g^{\beta_{o} \eta}\\)*
- After distribution, the original secret seeds are securely destroyed.

**Protocol (Proving)**

- Execute the program and get the assignment vector \\(v\\).
- Compute the linear-combinations of polynomials
  - \\(\ell(x) = \sum_{i=1}^{d} v_i \ell_{i}(x),\quad r(x) = \sum_{i=1}^{d} v_i r_{i}(x),\quad o(x) = \sum_{i=1}^{d} v_i o_{i}(x)\\)
- Compute the quotient polynomial:
  - \\(h(x) = \frac{\ell(x) r(x) - o(x)}{t(x)}\\)
- Evaluate each polynomial at \\(s\\):
  - \\(g^{\ell(s)} = \prod^{d}_{i=1} (g^{\ell_i(s)})^{v _i} ,\quad g^{r(s)} = \prod^{d} _{i=1} (g^{r_i(s)})^{v_i} ,\quad g^{o(s)} = \prod^{d} _{i=1} (g^{o _i(s)})^{v _i} \\)
- Evaluate each shifted polynomial at \\(s\\):
  - \\(g^{\alpha _{\ell} \ell(s)} = \prod^{d} _{i=1} (g^{\alpha _{\ell} \ell _i(s)})^{v _i} ,g^{\alpha _{r} r(s)} = \prod^{d} _{i=1} (g^{\alpha _{r} r _i(s)})^{v _i} ,g^{\alpha _{o} o(s)} = \prod^{d} _{i=1} (g^{\alpha _{o} o _i(s)})^{v _i} \\)
- Evaluate each consistency polynomial at \\(s\\):
  - \\(g^{z(s)} = \prod^{d}_{i=1} (g^{\beta _{\ell} \ell _{i}(s) + \beta _{r} r _{i}(s) + \beta _{o} o _{i}(s)})^{v _{i}}\\)
- Calculate \\(g^{h(s)}\\) from \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Proof**: 
  - \\((\\) \\(g^{\ell(s)}, g^{r(s)}, g^{o(s)}, g^{\alpha_{\ell} \ell(s)}, g^{\alpha_{r} r(s)}, g^{\alpha_{o} o(s)}, g^{h(s)},\\) \\(g^{z(s)}\\) \\()\\)

**Protocol (Verification)**

- Parse proof as \\((g^{\ell}, g^r, g^o, g^{\ell'}, g^{r'}, g^{o'}, g^{h}, g^{z})\\)
- Check polynomial restrictions
  - \\(e(g^{\ell}, g^{\alpha_{\ell}}) = e(g^{\ell'}, g),\quad e(g^{r}, g^{\alpha_{r}}) = e(g^{r'}, g),\quad e(g^{o}, g^{\alpha_{o}}) = e(g^{o'}, g)\\)
- Check variable consistency
  - *\\(e(g^{\ell}, g^{\beta _{\ell} \eta}) \cdot e(g^{r}, g^{\beta _{r} \eta}) \cdot e(g^{o}, g^{\beta _{o} \eta}) = e(g^{z}, g^{\eta})\\)*
- Validity check
  - \\(e(g^{\ell}, g^{r}) = e(g^t,g^h) \cdot e(g^o, g)\\)

**Implementation**

```rust
#[derive(Debug, Clone)]
pub struct ProofKey4 {
    g1_ell_i_vec: Vec<G1Point>,
    g1_r_i_vec: Vec<G1Point>,
    g2_r_i_vec: Vec<G2Point>,
    g1_o_i_vec: Vec<G1Point>,
    g1_alpha_ell_i_vec: Vec<G1Point>,
    g1_alpha_r_i_vec: Vec<G1Point>,
    g1_alpha_o_i_vec: Vec<G1Point>,
    g1_sj_vec: Vec<G1Point>,
    g1_checksum_vec: Vec<G1Point>,
}

#[derive(Debug, Clone)]
pub struct VerificationKey4 {
    g2_alpha_ell: G2Point,
    g2_alpha_r: G2Point,
    g2_alpha_o: G2Point,
    g2_beta_ell_eta: G2Point,
    g2_beta_r_eta: G2Point,
    g2_beta_o_eta: G2Point,
    g2_t_s: G2Point,
    g2_eta: G2Point,
}

#[derive(Debug, Clone)]
pub struct Proof4 {
    g1_ell: G1Point,
    g1_r: G1Point,
    g2_r: G2Point,
    g1_o: G1Point,
    g1_ell_prime: G1Point,
    g1_r_prime: G1Point,
    g1_o_prime: G1Point,
    g1_h: G1Point,
    g1_z: G1Point,
}

pub fn setup(g1: &G1Point, g2: &G2Point, qap: &QAP<FqOrder>) -> (ProofKey4, VerificationKey4) {
    let s = FqOrder::random_element(&[]);
    let alpha_ell = FqOrder::random_element(&[]);
    let alpha_r = FqOrder::random_element(&[]);
    let alpha_o = FqOrder::random_element(&[]);
    let beta_ell = FqOrder::random_element(&[]);
    let beta_r = FqOrder::random_element(&[]);
    let beta_o = FqOrder::random_element(&[]);
    let eta = FqOrder::random_element(&[]);

    let mut g1_checksum_vec = Vec::with_capacity(qap.d);

    for i in 0..qap.d {
        let ell_i_s = qap.ell_i_vec[i].eval(&s).sanitize();
        let r_i_s = qap.r_i_vec[i].eval(&s).sanitize();
        let o_i_s = qap.o_i_vec[i].eval(&s).sanitize();
        g1_checksum_vec.push(
            g1.mul_ref(
                ((beta_ell.mul_ref(&ell_i_s))
                    .add_ref(&beta_r.mul_ref(&r_i_s))
                    .add_ref(&beta_o.mul_ref(&o_i_s)))
                .get_value(),
            ),
        );
    }

    (
        ProofKey4 {
            g1_ell_i_vec: generate_challenge_vec(g1, &qap.ell_i_vec, &s),
            g1_r_i_vec: generate_challenge_vec(g1, &qap.r_i_vec, &s),
            g2_r_i_vec: generate_challenge_vec(g2, &qap.r_i_vec, &s),
            g1_o_i_vec: generate_challenge_vec(g1, &qap.o_i_vec, &s),
            g1_alpha_ell_i_vec: generate_alpha_challenge_vec(g1, &qap.ell_i_vec, &s, &alpha_ell),
            g1_alpha_r_i_vec: generate_alpha_challenge_vec(g1, &qap.r_i_vec, &s, &alpha_r),
            g1_alpha_o_i_vec: generate_alpha_challenge_vec(g1, &qap.o_i_vec, &s, &alpha_o),
            g1_sj_vec: generate_s_powers(g1, &s, qap.m),
            g1_checksum_vec: g1_checksum_vec,
        },
        VerificationKey4 {
            g2_alpha_ell: g2.mul_ref(alpha_ell.get_value()),
            g2_alpha_r: g2.mul_ref(alpha_r.get_value()),
            g2_alpha_o: g2.mul_ref(alpha_o.get_value()),
            g2_beta_ell_eta: g2.mul_ref(beta_ell.get_value()).mul_ref(eta.get_value()),
            g2_beta_r_eta: g2.mul_ref(beta_r.get_value()).mul_ref(eta.get_value()),
            g2_beta_o_eta: g2.mul_ref(beta_o.get_value()).mul_ref(eta.get_value()),
            g2_t_s: g2.mul_ref(qap.t.eval(&s).sanitize().get_value()),
            g2_eta: g2.mul_ref(eta.get_value()),
        },
    )
}

pub fn prove(assignment: &Vec<FqOrder>, proof_key: &ProofKey4, qap: &QAP<FqOrder>) -> Proof4 {
    Proof4 {
        g1_ell: accumulate_curve_points(&proof_key.g1_ell_i_vec, assignment),
        g1_r: accumulate_curve_points(&proof_key.g1_r_i_vec, assignment),
        g2_r: accumulate_curve_points(&proof_key.g2_r_i_vec, assignment),
        g1_o: accumulate_curve_points(&proof_key.g1_o_i_vec, assignment),
        g1_ell_prime: accumulate_curve_points(&proof_key.g1_alpha_ell_i_vec, assignment),
        g1_r_prime: accumulate_curve_points(&proof_key.g1_alpha_r_i_vec, assignment),
        g1_o_prime: accumulate_curve_points(&proof_key.g1_alpha_o_i_vec, assignment),
        g1_h: get_h(qap, assignment).eval_with_powers_on_curve(&proof_key.g1_sj_vec),
        g1_z: accumulate_curve_points(&proof_key.g1_checksum_vec, assignment),
    }
}

pub fn verify(
    g1: &G1Point,
    g2: &G2Point,
    proof: &Proof4,
    verification_key: &VerificationKey4,
) -> bool {
    let pairing1 = optimal_ate_pairing(&proof.g1_ell, &verification_key.g2_alpha_ell);
    let pairing2 = optimal_ate_pairing(&proof.g1_ell_prime, &g2);
    if pairing1 != pairing2 {
        return false;
    }

    let pairing3 = optimal_ate_pairing(&proof.g1_r, &verification_key.g2_alpha_r);
    let pairing4 = optimal_ate_pairing(&proof.g1_r_prime, &g2);
    if pairing3 != pairing4 {
        return false;
    }

    let pairing5 = optimal_ate_pairing(&proof.g1_o, &verification_key.g2_alpha_o);
    let pairing6 = optimal_ate_pairing(&proof.g1_o_prime, &g2);
    if pairing5 != pairing6 {
        return false;
    }

    let pairing7 = optimal_ate_pairing(&proof.g1_ell, &proof.g2_r);
    let pairing8 = optimal_ate_pairing(&proof.g1_h, &verification_key.g2_t_s);
    let pairing9 = optimal_ate_pairing(&proof.g1_o, &g2);

    if pairing7 != pairing8 * pairing9 {
        return false;
    }

    let pairing10 = optimal_ate_pairing(&proof.g1_ell, &verification_key.g2_beta_ell_eta);
    let pairing11 = optimal_ate_pairing(&proof.g1_r, &verification_key.g2_beta_r_eta);
    let pairing12 = optimal_ate_pairing(&proof.g1_o, &verification_key.g2_beta_o_eta);
    let pairing13 = optimal_ate_pairing(&proof.g1_z, &verification_key.g2_eta);

    pairing10 * pairing11 * pairing12 == pairing13
}
```

## Fifth Protocol: Pinocchio

The above protocol introduces four expensive pairing operations. To make it faster, Pinocchio protocl randomizes the generators.

**Protocol (Setup)**

- **Interpolated Polynomial:** Construct \\(\\{\ell_i, r_i, o_i\\}_{i\in[d]}\\) from \\(L\\), \\(R\\), and \\(O\\), respectively.
- **Target Polynomial:** \\(t(x) = (x-1)(x-2) \cdots (x-m)\\)
- **Secret Seed:** A trusted party generates the random value \\(s\\), \\(\alpha_{\ell}\\), \\(\alpha_r\\), \\(\alpha_o\\), \\(\beta\\), \\(\eta\\), *\\(\rho _{\ell}\\), and \\(\rho _{r}\\), and set \\(\rho _{o} = \rho _{\ell} \rho _{r}\\)*.
- **Randized Generators:** *\\(g _{\ell} = g^{\rho _{\ell}}\\), \\(g _{r} = g^{\rho _{r}}\\), and \\(g _{o} = g^{\rho _{o}}\\)*
- **Proof Key:** Provided to the prover
  - *\\(\\{g_{\ell}^{\ell_i(s)},g_{r}^{r_i(s)},g_{o}^{o_i(s)}\\}_{i\in[d]}\\)*
  - *\\(\\{g_{\ell}^{\alpha_{\ell} \ell_i(s)},g_{r}^{\alpha_{r} r_i(s)},g_{o}^{\alpha_{o} o_i(s)}\\}_{i\in[d]}\\)*
  - \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
  - *\\(\\{g _{\ell}^{\beta \ell _{i}(s)} \cdot g _{r}^{\beta r _{i}(s)} \cdot g _{o}^{\beta o _{i}(s)}\\} _{i \in [d]}\\)*
- **Verification Key:**
  - *\\(g _{o}^{t(s)}\\)* \\(, g^{\alpha _{\ell}}, g^{\alpha _{r}}, g^{\alpha _{o}}\\) *\\(,g^{\beta \eta}\\)*
- After distribution, the original secret seeds are securely destroyed.

**Protocol (Proving)**

- Execute the program and get the assignment vector \\(v\\).
- Compute the linear-combinations of polynomials
  - \\(\ell(x) = \sum_{i=1}^{d} v_i \ell_{i}(x),\quad r(x) = \sum_{i=1}^{d} v_i r_{i}(x),\quad o(x) = \sum_{i=1}^{d} v_i o_{i}(x)\\)
- Compute the quotient polynomial:
  - \\(h(x) = \frac{\ell(x) r(x) - o(x)}{t(x)}\\)
- Evaluate each polynomial at \\(s\\):
  - \\(g^{\ell(s)} = \prod^{d}_{i=1} (g _{\ell}^{\ell_i(s)})^{v _i} ,\quad g^{r(s)} = \prod^{d} _{i=1} (g _{r}^{r_i(s)})^{v_i} ,\quad g^{o(s)} = \prod^{d} _{i=1} (g _{o}^{o _i(s)})^{v _i} \\)
- Evaluate each shifted polynomial at \\(s\\):
  - \\(g^{\alpha _{\ell} \ell(s)} = \prod^{d} _{i=1} (g _{\ell}^{\alpha _{\ell} \ell _i(s)})^{v _i} ,g^{\alpha _{r} r(s)} = \prod^{d} _{i=1} (g _{r}^{\alpha _{r} r _i(s)})^{v _i} ,g^{\alpha _{o} o(s)} = \prod^{d} _{i=1} (g _{o}^{\alpha _{o} o _i(s)})^{v _i} \\)
- Evaluate each consistency polynomial at \\(s\\):
  - *\\(g^{z(s)} = \prod^{d}_{i=1} (g _{\ell}^{\beta \ell _{i}(s)} \cdot g _{r}^{\beta r _{i}(s)} \cdot g _{o}^{\beta o _{i}(s)})^{v _{i}}\\)*
- Calculate \\(g^{h(s)}\\) from \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Proof**: 
  - \\((\\) \\(g^{\ell(s)}, g^{r(s)}, g^{o(s)}, g^{\alpha_{\ell} \ell(s)}, g^{\alpha_{r} r(s)}, g^{\alpha_{o} o(s)}, g^{h(s)},\\) \\(g^{z(s)}\\) \\()\\)

**Protocol (Verification)**

- Parse proof as \\((g _{\ell}^{\ell}, g _{r}^r, g _{o}^o, g _{\ell}^{\ell'}, g _{r}^{r'}, g _{o}^{o'}, g^{h}, g^{z})\\)
- Check polynomial restrictions
  - *\\(e(g _{\ell}^{\ell}, g^{\alpha _{\ell}}) = e(g _{\ell}^{\ell'}, g)\\)*, *\\(e(g _{r}^{r}, g^{\alpha _{r}}) = e(g _{r}^{r'}, g)\\)*, *\\(e(g _{o}^{o}, g^{\alpha _{o}}) = e(g _{o}^{o'}, g)\\)*
- Check variable consistency
  - *\\(e(g _{\ell}^{\ell} \cdot g _{r}^{r} \cdot g _{o}^{o}) = e(g^{z}, g^{\eta})\\)*
- Validity check
  - *\\(e(g _{\ell}^{\ell}, g _{r}^{r}) = e(g _{o}^t,g^h) \cdot e(g _{o}^o, g)\\)*

This protocol reduces two pairing operations.