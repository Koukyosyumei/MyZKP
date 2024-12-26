# Proving Single Polynomial

Before dealing with all of \\(\ell(x)\\), \\(r(x)\\), and \\(o(x)\\) at once, we design a protocol that allows the Prover \\(\mathcal{A}\\) to convince the Verifier \\(\mathcal{B}\\) that \\(\mathcal{A}\\) knows a specific polynomial. Let's denote this polynomial of degree \\(n\\) with coefficients in a finite field as:

\begin{equation}
    P(x) = c_0 + c_1 x + c_2 x^{2} + \cdots c_n x^{n}
\end{equation}

Assume \\(P(x)\\) has \\(n\\) roots, \\(a_1, a_2, \ldots, a_n \in \mathbb{F}\\), such that \\(P(x) = (x - a_1)(x - a_2)\cdots(x - a_n)\\). The Verifier \\(\mathcal{B}\\) knows \\(m < n\\) roots of \\(P(x)\\), namely \\(a_1, a_2, \ldots, a_m\\). Let \\(T(x) = (x - a_1)(x - a_2)\cdots(x - a_m)\\). Note that the Prover also knows \\(T(x)\\).

The Prover's objective is to convince the Verifier that \\(\mathcal{A}\\) knows a polynomial \\(H(x) = \frac{P(x)}{T(x)}\\).

## Naive Approach

The simplest approach to prove that \\(\mathcal{A}\\) knows \\(H(x)\\) is as follows:

**Protocol:**

- \\(\mathcal{B}\\) sends all possible values in \\(\mathbb{F}\\) to \\(\mathcal{A}\\).
- \\(\mathcal{A}\\) computes and sends all possible outputs of \\(H(x)\\) and \\(P(x)\\).
- \\(\mathcal{B}\\) checks whether \\(H(a)T(a) = P(a)\\) holds for any \\(a\\) in \\(\mathbb{F}\\).

This protocol is highly inefficient, requiring \\(\mathcal{O}(|\mathbb{F}|)\\) evaluations and communications.

**Implementation:**

```rust
pub struct Prover<F: Field> {
    pub p: Polynomial<F>,
    pub t: Polynomial<F>,
    pub h: Polynomial<F>,
}

pub struct Verifier<F: Field> {
    pub t: Polynomial<F>,
    pub known_roots: Vec<F>,
}

impl<F: Field> Prover<F> {
    pub fn new(p: Polynomial<F>, t: Polynomial<F>) -> Self {
        let h = p.clone() / t.clone();
        Prover { p, t, h }
    }

    pub fn compute_all_values(&self, modulus: i128) -> (HashMap<F, F>, HashMap<F, F>) {
        let mut h_values = HashMap::new();
        let mut p_values = HashMap::new();

        for i in 0..modulus {
            let x = F::from_value(i);
            h_values.insert(x.clone(), self.h.eval(&x));
            p_values.insert(x.clone(), self.p.eval(&x));
        }

        (h_values, p_values)
    }
}

impl<F: Field> Verifier<F> {
    pub fn new(known_roots: Vec<F>) -> Self {
        let t = Polynomial::from_monomials(&known_roots);
        Verifier { t, known_roots }
    }

    pub fn verify(&self, h_values: &HashMap<F, F>, p_values: &HashMap<F, F>) -> bool {
        for (x, h_x) in h_values {
            let t_x = self.t.eval(x);
            let p_x = p_values.get(x).unwrap();
            if h_x.clone() * t_x != *p_x {
                return false;
            }
        }
        true
    }
}

pub fn naive_protocol<F: Field>(prover: &Prover<F>, verifier: &Verifier<F>, modulus: i128) -> bool {
    // Step 1: Verifier sends all possible values (implicitly done by Prover computing all values)

    // Step 2: Prover computes and sends all possible outputs
    let (h_values, p_values) = prover.compute_all_values(modulus);

    // Step 3: Verifier checks whether H(a)T(a) = P(a) holds for any a in F
    verifier.verify(&h_values, &p_values)
}
```

## \\(+\\) Schwartz-Zippel Lemma

Instead of evaluating polynomials at all values in \\(\mathbb{F}\\), we can leverage the Schwartz-Zippel Lemma: if \\(H(s) = \frac{P(s)}{T(s)}\\) or equivalently \\(H(s)T(s) = P(s)\\) for a random element \\(s\\), we can conclude that \\(H(x) = \frac{P(x)}{T(x)}\\) with high probability. Thus, the Prover \\(\mathcal{A}\\) only needs to send evaluations of \\(P(s)\\) and \\(H(s)\\) for a random input \\(s\\) received from \\(\mathcal{B}\\).

**Protocol:**
    
- *\\(\mathcal{B}\\) draws random \\(s\\) from \\(\mathbb{F}\\) and sends it to \\(\mathcal{A}\\).*
- *\\(\mathcal{A}\\) computes \\(h = H(s)\\) and \\(p = P(s)\\) and send them to \\(\mathcal{B}\\).*
- *\\(\mathcal{B}\\) checks whether \\(p = t h\\), where \\(t\\) denotes \\(T(s)\\).*

This protocol is efficient, requiring only a constant number of evaluations and communications.

**Implementation:**

```rust
pub struct Prover<F: Field> {
    pub p: Polynomial<F>,
    pub t: Polynomial<F>,
    pub h: Polynomial<F>,
}

pub struct Verifier<F: Field> {
    pub t: Polynomial<F>,
}

impl<F: Field> Prover<F> {
    pub fn new(p: Polynomial<F>, t: Polynomial<F>) -> Self {
        let h = p.clone() / t.clone();
        Prover { p, t, h }
    }

    pub fn compute_values(&self, s: &F) -> (F, F) {
        let h_s = self.h.eval(s);
        let p_s = self.p.eval(s);
        (h_s, p_s)
    }
}

impl<F: Field> Verifier<F> {
    pub fn new(t: Polynomial<F>) -> Self {
        Verifier { t }
    }

    pub fn generate_challenge(&self) -> F {
        F::random_element(&[])
    }

    pub fn verify(&self, s: &F, h: &F, p: &F) -> bool {
        let t_s = self.t.eval(s);
        h.clone() * t_s == *p
    }
}

pub fn schwartz_zippel_protocol<F: Field>(prover: &Prover<F>, verifier: &Verifier<F>) -> bool {
    // Step 1: Verifier generates a random challenge
    let s = verifier.generate_challenge();

    // Step 2: Prover computes and sends h and p
    let (h, p) = prover.compute_values(&s);

    // Step 3: Verifier checks whether p = t * h
    verifier.verify(&s, &h, &p)
}
```

**Vulnerability:**

However, it is vulnerable to a malicious prover who could send an arbitrary value \\(h'\\) and the corresponding \\(p'\\) such that \\(p' = h't\\).

```rust
// Simulating a malicious prover
pub struct MaliciousProver<F: Field> {
    t: Polynomial<F>,
}

impl<F: Field> MaliciousProver<F> {
    pub fn new(t: Polynomial<F>) -> Self {
        MaliciousProver { t }
    }

    pub fn compute_malicious_values(&self, s: &F) -> (F, F) {
        let h_prime = F::random_element(&[]);
        let t_s = self.t.eval(s);
        let p_prime = h_prime.clone() * t_s;
        (h_prime, p_prime)
    }
}

pub fn malicious_schwartz_zippel_protocol<F: Field>(
    prover: &MaliciousProver<F>,
    verifier: &Verifier<F>,
) -> bool {
    // Step 1: Verifier generates a random challenge
    let s = verifier.generate_challenge();

    // Step 2: Malicious Prover computes and sends h' and p'
    let (h_prime, p_prime) = prover.compute_malicious_values(&s);

    // Step 3: Verifier checks whether p' = t * h'
    verifier.verify(&s, &h_prime, &p_prime)
}

```

## \\(+\\) Discrete Logarithm Assumption

To address this vulnerability, the Verifier must hide the randomly chosen input \\(s\\) from the Prover. This can be achieved using the discrete logarithm assumption: it is computationally hard to determine \\(s\\) from \\(\alpha\\), where \\(\alpha = g^s \bmod p\\). Thus, it's safe for the Verifier to send \\(\alpha\\), as the Prover cannot easily derive \\(s\\) from it.

An interesting property of polynomial exponentiation is:

\begin{align}
    g^{P(x)} &= g^{c_0 + c_1 x + c_2 x^{2} + \cdots c_n x^{n}} = g^{c_0} (g^{x})^{c_1}  (g^{(x^2)})^{c_2} \cdots (g^{(x^n)})^{c_n}
\end{align}

Instead of sending \\(s\\), the Verifier can send \\(g\\) and \\(\alpha_{i} = g^{(s^i)}\\) for \\(i = 1, \cdots n\\). BE CAREFUL THAT **\\(g^{(s^i)} \neq (g^s)^i\\)**. The Prover can still evaluate \\(g^p = g^{P(s)}\\) using these powers of \\(g\\):

\begin{equation}
    g^{p} = g^{P(s)} = g^{c_0} \alpha_{1}^{c_1} (\alpha_{2})^{c_2} \cdots (\alpha_{n})^{c_n}
\end{equation}

Similarly, the Prover can evaluate \\(g^h = g^{H(s)}\\). The Verifier can then check \\(p = ht \iff g^p = (g^h)^t\\). 

**Protocol:**

- \\(\mathcal{B}\\) randomly draw \\(s\\) from \\(\mathbb{F}\\).
- *\\(\mathcal{B}\\) computes and sends \\(\\{\alpha, \alpha_2, ..., \alpha_{n}\\}\\), where \\(\alpha_i= g^{(s^{i})}\\).*
- *\\(\mathcal{A}\\) computes and sends \\(u = g^{p}\\) and \\(v = g^{h}\\).*
- *\\(\mathcal{B}\\) checks whether \\(u = v^{t}\\).*

This approach prevents the Prover from obtaining \\(s\\) or \\(t = T(s)\\), making it impossible to send fake \\((h', p')\\) such that \\(p' = h't\\).

```rust
pub struct Prover<F: Field> {
    pub p: Polynomial<F>,
    pub t: Polynomial<F>,
    pub h: Polynomial<F>,
}

pub struct Verifier<F: Field> {
    t: Polynomial<F>,
    s: F,
    g: F,
}

impl<F: Field> Prover<F> {
    pub fn new(p: Polynomial<F>, t: Polynomial<F>) -> Self {
        let h = p.clone() / t.clone();
        Prover { p, t, h }
    }

    pub fn compute_values(&self, alpha_powers: &[F]) -> (F, F) {
        let g_p = self.p.eval_with_powers(alpha_powers);
        let g_h = self.h.eval_with_powers(alpha_powers);
        (g_p, g_h)
    }
}

impl<F: Field> Verifier<F> {
    pub fn new(t: Polynomial<F>, generator: i128) -> Self {
        let mut rng = rand::thread_rng();
        let s = F::from_value(rng.gen_bigint_range(
            &BigInt::zero(),
            &BigInt::from(std::u64::MAX),
        ));
        let g = F::from_value(generator);
        Verifier { t, s, g }
    }

    pub fn generate_challenge(&self, max_degree: usize) -> Vec<F> {
        let mut alpha_powers = vec![];
        for i in 0..(max_degree + 1) {
            alpha_powers.push(
                self.g
                    .pow(self.s.clone().pow(i.to_bigint().unwrap()).get_value()),
            );
        }
        alpha_powers
    }

    pub fn verify(&self, u: &F, v: &F) -> bool {
        let t_s = self.t.eval(&self.s);
        u == &v.pow(t_s.get_value())
    }
}

pub fn discrete_log_protocol<F: Field>(prover: &Prover<F>, verifier: &Verifier<F>) -> bool {
    // Step 1 & 2: Verifier generates a challenge
    let max_degree = prover.p.degree();
    let alpha_powers = verifier.generate_challenge(max_degree as usize);

    // Step 3: Prover computes and sends u = g^p and v = g^h
    let (u, v) = prover.compute_values(&alpha_powers);

    // Step 4: Verifier checks whether u = v^t
    verifier.verify(&u, &v)
}
```

**Vulnerability:**

However, this protocol still has a flaw. Since the Prover can compute \\(g^t = \alpha_{c_1}(\alpha_2)^{c_2}\cdots(\alpha_m)^{c_m}\\), they could send fake values \\(((g^{t})^{z}, g^{z})\\) instead of \\((g^p, g^h)\\) for an arbitrary value \\(z\\). The verifier's check would still pass, and they could not detect this deception.

```rust
// Simulating a malicious prover
pub struct MaliciousProver<F: Field> {
    t: Polynomial<F>,
}

impl<F: Field> MaliciousProver<F> {
    pub fn new(t: Polynomial<F>) -> Self {
        MaliciousProver { t }
    }

    pub fn compute_malicious_values(&self, alpha_powers: &[F]) -> (F, F) {
        let g_t = self.t.eval_with_powers(alpha_powers);
        let z = F::random_element(&[]);
        let g = &alpha_powers[0];
        let fake_v = g.pow(z.get_value());
        let fake_u = g_t.pow(z.get_value());
        (fake_u, fake_v)
    }
}

pub fn malicious_discrete_log_protocol<F: Field>(
    prover: &MaliciousProver<F>,
    verifier: &Verifier<F>,
) -> bool {
    // Step 1 & 2: Verifier generates a challenge
    let max_degree = prover.t.degree() as usize;
    let alpha_powers = verifier.generate_challenge(max_degree as usize);

    // Step 3: Malicious Prover computes and sends fake u and v
    let (fake_u, fake_v) = prover.compute_malicious_values(&alpha_powers);

    // Step 4: Verifier checks whether u = v^t (which will pass for the fake values)
    verifier.verify(&fake_u, &fake_v)
}
```

## \\(+\\) Knowledge of Exponent Assumption

To address the vulnerability where the verifier \\(\mathcal{B}\\) cannot distinguish if \\(v (= g^h)\\) from the prover is a power of \\(\alpha_i = g^{(s^i)}\\), we can employ the Knowledge of Exponent Assumption. This approach involves the following steps:

- \\(\mathcal{B}\\) sends both \\(\alpha_i\\) and \\(\alpha'_i = \alpha_i^r\\) for a new random value \\(r\\).
- The prover returns \\(a = (\alpha_i)^{c_i}\\) and \\(a' = (\alpha'\_i)^{c_i}\\) for \\(i = 1, ..., n\\).
- \\(\mathcal{B}\\) can conclude that \\(a\\) is a power of \\(\alpha_i\\) if \\(a^r = a'\\).


Based on this assumption, we can design an improved protocol:

**Protocol:**

- \\(\mathcal{B}\\) randomly selects \\(s\\) and *\\(r\\)* from field \\(\mathbb{F}\\).
- \\(\mathcal{B}\\) computes and sends \\(\\{\alpha_1, \alpha_2, ..., \alpha_{n}\\}\\) *and \\(\\{\alpha'\_1, \alpha'\_2, ..., \alpha'\_{n}\\}\\), where \\(\alpha_i = g^{(s^i)}\\) and \\(\alpha' = \alpha_{r} = g^{(s^{i})r}\\).*
- \\(\mathcal{A}\\) computes and sends \\(u = g^{p}\\), \\(v = g^{h}\\), *and \\(w = g^{p'}\\), where \\(g^{p'} = g^{rP(s)}\\).*
- *\\(\mathcal{B}\\) checks whether \\(u^{r} = w\\).*
- \\(\mathcal{B}\\) checks whether \\(u = v^{t}\\).

The prover can compute \\(g^{p'} = g^{rP(s)} = \alpha'^{c_1} (\alpha'^{2})^{c_2} \cdots (\alpha'^{n})^{c_n}\\) using powers of \\(\alpha'\\). This protocol now satisfies the properties of a SNARK: completeness, soundness, and efficiency.

**Implementation:**

```rust
pub struct Prover<F: Field> {
    pub p: Polynomial<F>,
    pub t: Polynomial<F>,
    pub h: Polynomial<F>,
}

pub struct Verifier<F: Field> {
    t: Polynomial<F>,
    s: F,
    r: F,
    g: F,
}

impl<F: Field> Prover<F> {
    pub fn new(p: Polynomial<F>, t: Polynomial<F>) -> Self {
        let h = p.clone() / t.clone();
        Prover { p, t, h }
    }

    pub fn compute_values(&self, alpha_powers: &[F], alpha_prime_powers: &[F]) -> (F, F, F) {
        let g_p = self.p.eval_with_powers(alpha_powers);
        let g_h = self.h.eval_with_powers(alpha_powers);
        let g_p_prime = self.p.eval_with_powers(alpha_prime_powers);
        (g_p, g_h, g_p_prime)
    }
}

impl<F: Field> Verifier<F> {
    pub fn new(t: Polynomial<F>, generator: i128) -> Self {
        let mut rng = rand::thread_rng();
        let s = F::from_value(rng.gen_bigint_range(
            &BigInt::zero(),
            &BigInt::from(std::u64::MAX),
        ));
        let r = F::from_value(rng.gen_bigint_range(
            &BigInt::zero(),
            &BigInt::from(std::u64::MAX),
        ));
        let g = F::from_value(generator);
        Verifier { t, s, r, g }
    }

    pub fn generate_challenge(&self, max_degree: usize) -> (Vec<F>, Vec<F>) {
        let mut alpha_powers = vec![];
        let mut alpha_prime_powers = vec![];

        for i in 0..(max_degree + 1) {
            alpha_powers.push(
                self.g
                    .pow(self.s.clone().pow(i.to_bigint().unwrap()).get_value()),
            );
            alpha_prime_powers.push(alpha_powers.last().unwrap().pow(self.r.get_value()));
        }

        (alpha_powers, alpha_prime_powers)
    }

    pub fn verify(&self, u: &F, v: &F, w: &F) -> bool {
        let t_s = self.t.eval(&self.s);
        let u_r = u.pow(self.r.clone().get_value());

        // Check 1: u^r = w
        let check1 = u_r == *w;

        // Check 2: u = v^t
        let check2 = *u == v.pow(t_s.get_value());

        check1 && check2
    }
}

pub fn knowledge_of_exponent_protocol<F: Field>(
    prover: &Prover<F>,
    verifier: &Verifier<F>,
) -> bool {
    // Step 1 & 2: Verifier generates a challenge
    let max_degree = std::cmp::max(prover.p.degree(), prover.h.degree()) as usize;
    let (alpha_powers, alpha_prime_powers) = verifier.generate_challenge(max_degree + 1);

    // Step 3: Prover computes and sends u = g^p, v = g^h, and w = g^p'
    let (u, v, w) = prover.compute_values(&alpha_powers, &alpha_prime_powers);

    // Step 4 & 5: Verifier checks whether u^r = w and u = v^t
    verifier.verify(&u, &v, &w)
}
```

## \\(+\\) Zero Knowledge

To transform the above protocol into a zk-SNARK, we need to ensure that the verifier cannot learn anything about \\(P(x)\\) from the prover's information. This is achieved by having the prover obfuscate all information with a random secret value \\(\delta\\):

**Protocol:**

- \\(\mathcal{B}\\) randomly selects \\(s\\) and \\(r\\) from field \\(\mathbb{F}\\).
- \\(\mathcal{B}\\) computes and sends \\(\\{\alpha_1, \alpha_2, ..., \alpha_{n}\\}\\) and \\(\\{\alpha\_1', \alpha'\_2, ..., \alpha'\_{n}\\}\\), where \\(\alpha_i = g^{(s^{i})}\\) and \\(\alpha_i' = \alpha_i^{r} = g^{(s^{i})r}\\).
- *\\(\mathcal{A}\\) randomly selects \\(\delta\\) from field \\(\mathbb{F}\\).*
- *\\(\mathcal{A}\\) computes and sends \\(u' = (g^{p})^{\delta}\\), \\(v' = (g^{h})^{\delta}\\), and \\(w' = (g^{p'})^{\delta}\\).*
- \\(\mathcal{B}\\) checks whether \\(u'^{r} = w'\\).
- \\(\mathcal{B}\\) checks whether \\(u' = v'^{t}\\).

By introducing the random value \\(\delta\\), the verifier can no longer learn anything about \\(p\\), \\(h\\), or \\(w\\), thus achieving zero knowledge.


**Implementation:**

```rust
pub struct Prover<F: Field> {
    pub p: Polynomial<F>,
    pub t: Polynomial<F>,
    pub h: Polynomial<F>,
}

pub struct Verifier<F: Field> {
    t: Polynomial<F>,
    s: F,
    r: F,
    g: F,
}

impl<F: Field> Prover<F> {
    pub fn new(p: Polynomial<F>, t: Polynomial<F>) -> Self {
        let h = p.clone() / t.clone();
        Prover { p, t, h }
    }

    pub fn compute_values(&self, alpha_powers: &[F], alpha_prime_powers: &[F]) -> (F, F, F) {
        let mut rng = rand::thread_rng();
        let delta = F::from_value(rng.gen_bigint_range(
            &BigInt::zero(),
            &BigInt::from(std::u64::MAX),
        ));

        let g_p = self.p.eval_with_powers(alpha_powers);
        let g_h = self.h.eval_with_powers(alpha_powers);
        let g_p_prime = self.p.eval_with_powers(alpha_prime_powers);

        let u_prime = g_p.pow(delta.get_value());
        let v_prime = g_h.pow(delta.get_value());
        let w_prime = g_p_prime.pow(delta.get_value());

        (u_prime, v_prime, w_prime)
    }
}

impl<F: Field> Verifier<F> {
    pub fn new(t: Polynomial<F>, generator: i128) -> Self {
        let mut rng = rand::thread_rng();
        let s = F::from_value(rng.gen_bigint_range(
            &BigInt::zero(),
            &BigInt::from(std::u64::MAX),
        ));
        let r = F::from_value(rng.gen_bigint_range(
            &BigInt::zero(),
            &BigInt::from(std::u64::MAX),
        ));
        let g = F::from_value(generator);
        Verifier { t, s, r, g }
    }

    pub fn generate_challenge(&self, max_degree: usize) -> (Vec<F>, Vec<F>) {
        let mut alpha_powers = vec![];
        let mut alpha_prime_powers = vec![];

        for i in 0..(max_degree + 1) {
            alpha_powers.push(
                self.g
                    .pow(self.s.clone().pow(i.to_bigint().unwrap()).get_value()),
            );
            alpha_prime_powers.push(alpha_powers.last().unwrap().pow(self.r.get_value()));
        }

        (alpha_powers, alpha_prime_powers)
    }

    pub fn verify(&self, u_prime: &F, v_prime: &F, w_prime: &F) -> bool {
        let t_s = self.t.eval(&self.s);
        let u_prime_r = u_prime.pow(self.r.clone().get_value());

        // Check 1: u'^r = w'
        let check1 = u_prime_r == *w_prime;

        // Check 2: u' = v'^t
        let check2 = *u_prime == v_prime.pow(t_s.get_value());

        check1 && check2
    }
}

pub fn zk_protocol<F: Field>(prover: &Prover<F>, verifier: &Verifier<F>) -> bool {
    // Step 1 & 2: Verifier generates a challenge
    let max_degree = std::cmp::max(prover.p.degree(), prover.h.degree()) as usize;
    let (alpha_powers, alpha_prime_powers) = verifier.generate_challenge(max_degree + 1);

    // Step 3 & 4: Prover computes and sends u' = (g^p)^δ, v' = (g^h)^δ, and w' = (g^p')^δ
    let (u_prime, v_prime, w_prime) = prover.compute_values(&alpha_powers, &alpha_prime_powers);

    // Step 5 & 6: Verifier checks whether u'^r = w' and u' = v'^t
    verifier.verify(&u_prime, &v_prime, &w_prime)
}
```

## \\(+\\) Non-interactivity

The previously described protocol requires each verifier to generate unique random values, which becomes inefficient when a prover needs to demonstrate knowledge to multiple verifiers. To address this, we aim to eliminate the interaction between the prover and verifier. One effective solution is the use of a trusted setup.

**Protocol (Trusted Setup):**

- ***Secret Seed:** A trusted third party generates the random values \\(s\\) and \\(r\\)*
- ***Proof Key:** Provided to the prover*
    - *\\(\\{\alpha_1, \alpha_2, ..., \alpha_{n}\\}\\), where \\(\alpha_{i} = g^{(s^i)}\\)*
    - *\\(\\{\alpha'\_1, \alpha'\_2, ..., \alpha'\_{n}\\}\\), where \\(\alpha_i' = g^{(s^{i})r}\\)*
- ***Verification Key:** Distributed to verifiers*
    - *\\(g^{r}\\)*
    - *\\(g^{T(s)}\\)* 
- *After distribution, the original \\(s\\) and \\(r\\) values are securely destroyed.*


Then, the non-interactive protocol consists of two main parts: proof generation and verification.

**Protocol (Proof):**
    
- *\\(\mathcal{A}\\) receives the proof key*
- *\\(\mathcal{A}\\) randomly selects \\(\delta\\) from field \\(\mathbb{F}\\).*
- *\\(\mathcal{A}\\) broadcast the proof \\(\pi = (u' = (g^{p})^{\delta}, v' = (g^{h})^{\delta}, w' = (g^{p'})^{\delta})\\)*
    

Since \\(r\\) is not shared and already destroyed, the verifier \\(\mathcal{B}\\) cannot calculate \\(u'^{r}\\) to check \\(u'^{r} = w'\\). Instead, the verifier can use a paring with bilinear mapping; \\(u'^{r} = w'\\) is equivalent to \\(e(u' = (g^{p})^{\delta}, g^{r}) = e(w'=(g^{p'})^{\delta}, g)\\).

**Protocol (Verification):**
    
- *\\(\mathcal{B}\\) receives the verification key.*
- *\\(\mathcal{B}\\) receives the proof \\(\pi\\).*
- *\\(\mathcal{B}\\) checks whether \\(e(u', g^{r}) = e(w', g)\\).*
- *\\(\mathcal{B}\\) checks whether \\(e(u', g) = e (v', g^{T(s)})\\).*

**Implementation:**

```rust
pub struct ProofKey {
    alpha: Vec<G1Point>,
    alpha_prime: Vec<G1Point>,
}

pub struct VerificationKey {
    g_r: G2Point,
    g_t_s: G2Point,
}

pub struct Proof {
    u_prime: G1Point,
    v_prime: G1Point,
    w_prime: G1Point,
}

pub struct TrustedSetup {
    proof_key: ProofKey,
    verification_key: VerificationKey,
}

impl TrustedSetup {
    pub fn generate(g1: &G1Point, g2: &G2Point, t: &Polynomial<Fq>, n: usize) -> Self {
        let mut rng = rand::thread_rng();
        let s = Fq::from_value(rng.gen_bigint_range(&BigInt::zero(), &BigInt::from(std::u32::MAX)));
        let r = Fq::from_value(rng.gen_bigint_range(&BigInt::zero(), &BigInt::from(std::u32::MAX)));

        let mut alpha = Vec::with_capacity(n);
        let mut alpha_prime = Vec::with_capacity(n);

        let mut s_power = Fq::one();
        for _ in 0..1 + n {
            alpha.push(g1.clone() * s_power.clone().get_value());
            alpha_prime.push(g1.clone() * (s_power.clone() * r.clone()).get_value());
            s_power = s_power * s.clone();
        }

        let g_r = g2.clone() * r.clone().get_value();
        let g_t_s = g2.clone() * t.eval(&s).get_value();

        TrustedSetup {
            proof_key: ProofKey { alpha, alpha_prime },
            verification_key: VerificationKey { g_r, g_t_s },
        }
    }
}

pub struct Prover {
    pub p: Polynomial<Fq>,
    pub h: Polynomial<Fq>,
}

impl Prover {
    pub fn new(p: Polynomial<Fq>, t: Polynomial<Fq>) -> Self {
        let h = p.clone() / t;
        Prover { p, h }
    }

    pub fn generate_proof(&self, proof_key: &ProofKey) -> Proof {
        let mut rng = rand::thread_rng();
        let delta =
            Fq::from_value(rng.gen_bigint_range(&BigInt::zero(), &BigInt::from(std::u32::MAX)));

        let g_p = self.p.eval_with_powers_on_curve(&proof_key.alpha);
        let g_h = self.h.eval_with_powers_on_curve(&proof_key.alpha);
        let g_p_prime = self.p.eval_with_powers_on_curve(&proof_key.alpha_prime);

        Proof {
            u_prime: g_p * delta.get_value(),
            v_prime: g_h * delta.get_value(),
            w_prime: g_p_prime * delta.get_value(),
        }
    }
}

pub struct Verifier {
    pub g1: G1Point,
    pub g2: G2Point,
}

impl Verifier {
    pub fn new(g1: G1Point, g2: G2Point) -> Self {
        Verifier { g1, g2 }
    }

    pub fn verify(&self, proof: &Proof, vk: &VerificationKey) -> bool {
        // Check e(u', g^r) = e(w', g)
        let pairing1 = optimal_ate_pairing(&proof.u_prime, &vk.g_r);
        let pairing2 = optimal_ate_pairing(&proof.w_prime, &self.g2);
        let check1 = pairing1 == pairing2;

        // Check e(u', g^t) = e(v', g)
        let pairing3 = optimal_ate_pairing(&proof.u_prime, &self.g2);
        let pairing4 = optimal_ate_pairing(&proof.v_prime, &vk.g_t_s);
        let check2 = pairing3 == pairing4;

        check1 && check2
    }
}

pub fn non_interactive_zkp_protocol(
    prover: &Prover,
    verifier: &Verifier,
    setup: &TrustedSetup,
) -> bool {
    let proof = prover.generate_proof(&setup.proof_key);
    verifier.verify(&proof, &setup.verification_key)
}
```
    

