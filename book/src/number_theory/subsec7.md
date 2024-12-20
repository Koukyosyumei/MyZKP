# Assumptions

### Assumption: Discrete Logarithm Problem 

---

*Let \\(G\\) be a finite cyclic group of order \\(n\\), with \\(\gamma\\) as its generator and \\(1\\) as the identity element. For any element \\(\alpha \in G\\), there is currently no known efficient (polynomial-time) algorithm to compute the smallest non-negative integer \\(x\\) such that \\(\alpha = \gamma^{x}\\).*

---

The Discrete Logarithm Problem can be thought of as a one-way function. It's easy to compute \\(g^{x}\\) given \\(g\\) and \\(x\\), but it's computationally difficult to find \\(x\\) given \\(g\\) and \\(g^{x}\\).

### Assumption: Elliptic Curve Discrete Logarithm Problem

---

*Let \\(E\\) be an elliptic curve defined over a finite field \\(\mathbb{F}_q\\), where \\(q\\) is a prime power. Let \\(P\\) be a point on \\(E\\) of point order \\(n\\), and let \\(\langle P \rangle\\) be the cyclic subgroup of \\(E\\) generated by \\(P\\). For any element \\(Q \in \langle P \rangle\\), there is currently no known efficient (polynomial-time) algorithm to compute the unique integer \\(k\\), \\(0 \leq k < n\\), such that \\(Q = kP\\).*

---

This assumption is an elliptic curve version of the Discrete Logarithm Problem.

### Assumption: Knowledge of Exponent Assumption

---

*Let \\(G\\) be a cyclic group of prime order \\(q\\) with generator \\(g \in G\\). For any probabilistic polynomial-time algorithm \\(\mathcal{A}\\) that outputs:*

\begin{equation}
\mathcal{A}(g, g^x) = (h, h') \quad s.t. \quad h' = h^x
\end{equation}
*, there exists an efficient extractor \\(\mathcal{E}\\) such that:*
\begin{equation}
\mathcal{E}(\mathcal{A}, g, g^x) = y \quad s.t. \quad h = g^y
\end{equation}

---

This assumption states that if \\(\mathcal{A}\\) can compute the pair \\((g^y, g^{xy})\\) from \\((g, g^x)\\), then \\(\mathcal{A}\\) must "know" the value \\(y\\), in the sense that \\(\mathcal{E}\\) can extract \\(y\\) from \\(\mathcal{A}\\)'s internal state.
The Knowledge of Exponent Assumption is useful for constructing verifiable exponential calculation algorithms. Consider a scenario where Alice has a secret value \\(a\\), and Bob has a secret value \\(b\\). Bob wants to obtain \\(g^{ab}\\). This can be achieved through the following protocol:

**Verifiable Exponential Calculation Algorithm**

1. Bob sends \\((g, g'=g^{b})\\) to Alice
2. Alice sends \\((h=g^{a}, h'=g'^{a})\\) to Bob
3. Bob checks \\(h^{b} = h'\\).

Thanks to the Discrete Logarithm Assumption and the Knowledge of Exponent Assumption, the following properties hold:

- Bob cannot derive \\(a\\) from \\((h, h')\\).
- Alice cannot derive \\(b\\) from \\((g, g')\\).
- Alice cannot generate \\((t, t')\\) such that \\(t \neq h\\) and \\(t^{b} = t'\\).  
- If \\(h^{b} = h'\\), Bob can conclude that \\(h\\) is the power of \\(g\\).