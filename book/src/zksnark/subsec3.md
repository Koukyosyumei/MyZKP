# Proving Single Polynomial

Before dealing with all of \\(\ell(x)\\), \\(r(x)\\), and \\(o(x)\\) at once, we design a protocol that allows the Prover \\(\mathcal{A}\\) to convince the Verifier \\(\mathcal{B}\\) that \\(\mathcal{A}\\) knows a specific polynomial. Let's denote this polynomial of degree \\(n\\) with coefficients in a finite field as:

\begin{equation}
    P(x) = c_0 + c_1 x + c_2 x^{2} + \cdots c_n x^{n}
\end{equation}

Assume \\(P(x)\\) has \\(n\\) roots, \\(a_1, a_2, \ldots, a_n \in \mathbb{F}\\), such that \\(P(x) = (x - a_1)(x - a_2)\cdots(x - a_n)\\). The Verifier \\(\mathcal{B}\\) knows \\(d < n\\) roots of \\(P(x)\\), namely \\(a_1, a_2, \ldots, a_d\\). Let \\(T(x) = (x - a_1)(x - a_2)\cdots(x - a_d)\\). Note that the Prover also knows \\(T(x)\\).

The Prover's objective is to convince the Verifier that \\(\mathcal{A}\\) knows a polynomial \\(H(x) = \frac{P(x)}{T(x)}\\).

### Naive Approach

The simplest approach to prove that \\(\mathcal{A}\\) knows \\(H(x)\\) is as follows:

*Naive Approach*

- *\\(\mathcal{B}\\) sends all possible values in \\(\mathbb{F}\\) to \\(\mathcal{A}\\).*
- *\\(\mathcal{A}\\) computes and sends all possible outputs of \\(H(x)\\) and \\(P(x)\\).*
- *\\(\mathcal{B}\\) checks whether \\(H(a)T(a) = P(a)\\) holds for any \\(a\\) in \\(\mathbb{F}\\).*

This protocol is highly inefficient, requiring \\(\mathcal{O}(|\mathbb{F}|)\\) evaluations and communications.

### \\(+\\) Schwartz-Zippel Lemma

Instead of evaluating polynomials at all values in \\(\mathbb{F}\\), we can leverage the Schwartz-Zippel Lemma: if \\(H(s) = \frac{P(s)}{T(s)}\\) or equivalently \\(H(s)T(s) = P(s)\\) for a random element \\(s\\), we can conclude that \\(H(x) = \frac{P(x)}{T(x)}\\) with high probability. Thus, the Prover \\(\mathcal{A}\\) only needs to send evaluations of \\(P(s)\\) and \\(H(s)\\) for a random input \\(s\\) received from \\(\mathcal{B}\\).

*\\(+\\) Schwartz-Zippel Lemma*
    
- *\\(\mathcal{B}\\) draws random \\(s\\) from \\(\mathbb{F}\\) and sends it to \\(\mathcal{A}\\).*
- *\\(\mathcal{A}\\) computes \\(h = H(s)\\) and \\(p = P(s)\\) and send them to \\(\mathcal{B}\\).*
- *\\(\mathcal{B}\\) checks whether \\(p = t h\\), where \\(t\\) denotes \\(T(s)\\).*

This protocol is efficient, requiring only a constant number of evaluations and communications. However, it is vulnerable to a malicious prover who could send an arbitrary value \\(h'\\) and the corresponding \\(p'\\) such that \\(p' = h't\\).

### \\(+\\) Discrete Logarithm Assumption

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

*\\(+\\) Discrete Logarithm Assumption*

- *\\(\mathcal{B}\\) randomly draw \\(s\\) from \\(\mathbb{F}\\).*
- *\\(\mathcal{B}\\) computes and sends \\(\\{\alpha, \alpha_2, ..., \alpha_{n}\\}\\), where \\(\alpha_i= g^{(s^{i})}\\).*
- *\\(\mathcal{A}\\) computes and sends \\(u = g^{p}\\) and \\(v = g^{h}\\).*
- *\\(\mathcal{B}\\) checks whether \\(u = v^{t}\\).*

This approach prevents the Prover from obtaining \\(s\\) or \\(t = T(s)\\), making it impossible to send fake \\((h', p')\\) such that \\(p' = h't\\).

However, this protocol still has a flaw. Since the Prover can compute \\(g^t = \alpha_{c_1}(\alpha_2)^{c_2}\cdots(\alpha_d)^{c_d}\\), they could send fake values \\(((g^{t})^{z}, g^{z})\\) instead of \\((g^p, g^h)\\) for an arbitrary value \\(z\\). The verifier's check would still pass, and they could not detect this deception.

### \\(+\\) Knowledge of Exponent Assumption

To address the vulnerability where the verifier \\(\mathcal{B}\\) cannot distinguish if \\(v (= g^h)\\) from the prover is a power of \\(\alpha_i = g^{(s^i)}\\), we can employ the Knowledge of Exponent Assumption. This approach involves the following steps:

- \\(\mathcal{B}\\) sends both \\(\alpha_i\\) and \\(\alpha'_i = \alpha_i^r\\) for a new random value \\(r\\).
- The prover returns \\(a = (\alpha_i)^{c_i}\\) and \\(a' = (\alpha'\_i)^{c_i}\\) for \\(i = 1, ..., n\\).
- \\(\mathcal{B}\\) can conclude that \\(a\\) is a power of \\(\alpha_i\\) if \\(a^r = a'\\).


Based on this assumption, we can design an improved protocol:

*\\(+\\) Knowledge of Exponent Assumption*

- *\\(\mathcal{B}\\) randomly selects \\(s\\) and \\(r\\) from field \\(\mathbb{F}\\).*
- *\\(\mathcal{B}\\) computes and sends \\(\{\alpha_1, \alpha_2, ..., \alpha_{n}\}\\) and \\(\{\alpha'\_1, \alpha'\_2, ..., \alpha'\_{n}\}\\), where \\(\alpha_i = g^{(s^i)}\\) and \\(\alpha' = \alpha_{r} = g^{(s^{i})r}\\).*
- *\\(\mathcal{A}\\) computes and sends \\(u = g^{p}\\), \\(v = g^{h}\\), and \\(w = g^{p'}\\), where \\(g^{p'} = g^{P(sr)}\\).*
- *\\(\mathcal{B}\\) checks whether \\(u^{r} = w\\).*
- *\\(\mathcal{B}\\) checks whether \\(u = v^{t}\\).*


The prover can compute \\(g^{p'} = g^{P(sr)} = \alpha'^{c_1} (\alpha'^{2})^{c_2} \cdots (\alpha'^{n})^{c_n}\\) using powers of \\(\alpha'\\). This protocol now satisfies the properties of a SNARK: completeness, soundness, and efficiency.

### \\(+\\) Zero Knowledge

To transform the above protocol into a zk-SNARK, we need to ensure that the verifier cannot learn anything about \\(P(x)\\) from the prover's information. This is achieved by having the prover obfuscate all information with a random secret value \\(\delta\\):

*\\(+\\) Zero Knowledge*

- *\\(\mathcal{B}\\) randomly selects \\(s\\) and \\(r\\) from field \\(\mathbb{F}\\).*
- *\\(\mathcal{B}\\) computes and sends \\(\\{\alpha_1, \alpha_2, ..., \alpha_{n}\\}\\) and \\(\\{\alpha\_1', \alpha'\_2, ..., \alpha'\_{n}\\}\\), where \\(\alpha_i = g^{(s^{i})}\\) and \\(\alpha_i' = \alpha_i^{r} = g^{(s^{i})r}\\).*
- *\\(\mathcal{A}\\) randomly selects \\(\delta\\) from field \\(\mathbb{F}\\).*
- *\\(\mathcal{A}\\) computes and sends \\(u' = (g^{p})^{\delta}\\), \\(v' = (g^{h})^{\delta}\\), and \\(w' = (g^{p'})^{\delta}\\).*
- *\\(\mathcal{B}\\) checks whether \\(u'^{r} = w'\\).*
- *\\(\mathcal{B}\\) checks whether \\(u' = v'^{t}\\).*



By introducing the random value \\(\delta\\), the verifier can no longer learn anything about \\(p\\), \\(h\\), or \\(w\\), thus achieving zero knowledge.

### \\(+\\) Non-interactivity

The previously described protocol requires each verifier to generate unique random values, which becomes inefficient when a prover needs to demonstrate knowledge to multiple verifiers. To address this, we aim to eliminate the interaction between the prover and verifier. One effective solution is the use of a trusted setup.

*\\(+\\) Non-Interactive: Trusted Setup*

- ***Secret Seed:** A trusted third party generates the random values \\(s\\) and \\(r\\)*
- ***Proof Key:** Provided to the prover*
    - *\\(\\{\alpha_1, \alpha_2, ..., \alpha_{n}\\}\\), where \\(\alpha_{i} = g^{(s^i)}\\)*
    - *\\(\\{\alpha'\_1, \alpha'\_2, ..., \alpha'\_{n}\\}\\), where \\(\alpha_i' = g^{(s^{i})r}\\)*
- ***Verification Key:** Distributed to verifiers*
    - *\\(g^{r}\\)*
    - *\\(g^{T(s)}\\)* 
- *After distribution, the original \\(s\\) and \\(r\\) values are securely destroyed.*


Then, the non-interactive protocol consists of two main parts: proof generation and verification.

*\\(+\\) Non-Interactive: Proof*
    
- *\\(\mathcal{A}\\) receives the proof key*
- *\\(\mathcal{A}\\) randomly selects \\(\delta\\) from field \\(\mathbb{F}\\).*
- *\\(\mathcal{A}\\) broadcast the proof \\(\pi = (u' = (g^{p})^{\delta}, v' = (g^{h})^{\delta}, w' = (g^{p'})^{\delta})\\)*
    


Since \\(r\\) is not shared and already destroyed, the verifier \\(\mathcal{B}\\) cannot calculate \\(u'^{r}\\) to check \\(u'^{r} = w'\\). Instead, the verifier can use a paring with bilinear mapping; \\(u'^{r} = w'\\) is equivalent to \\(e(u' = (g^{p})^{\delta}, g^{r}) = e(w'=(g^{p'})^{\delta}, g)\\).

*\\(+\\) Non-Interactive: Verification*
    
- *\\(\mathcal{B}\\) receives the verification key.*
- *\\(\mathcal{B}\\) receives the proof \\(\pi\\).*
- *\\(\mathcal{B}\\) checks whether \\(e(u', g^{r}) = e(w', g)\\).*
- *\\(\mathcal{B}\\) checks whether \\(u' = v'^{t}\\).*
    

