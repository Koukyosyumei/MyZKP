# FRI

The Fast Reed-Solomon IOP of Proximity (FRI) is a powerful interactive proof system that allows a verifier to efficiently check if a claimed polynomial has a degree bounded by some value \\(d\\), without the prover needing to reveal the entire polynomial.

## Core Concept

At its heart, FRI employs a clever divide-and-conquer strategy. The prover commits to the evaluation of the polynomial over a large domain. The protocol then iteratively "folds" the polynomial into a new polynomial of roughly half the degree. This process continues until the degree is so small that the verifier can directly check the degree bound by querying the remaining evaluations. The logarithmic communication complexity, scaling as \\(\mathcal{O}(\log{d})\\), arises from this halving of the degree in each round.

## Mathematical Foundation

Consider a polynomial \\(f(X) = \sum^{d}_{i=0} c_i X^{i}\\), where \\(d\\) is the maximum degree. This polynomial can be divided into two parts:

\\[f(X) = f_E(X^2) + X \cdot f_O(X^2)\\]

, where 

- \\(f\_E(X^2) = \frac{f(X) + f(-X)}{2} = \sum\_{i=0}^{\frac{d + 1}{2} - 1} c\_{2i} X^{2i} \\)
- \\(f\_O(X^2) = \frac{f(X) - f(-X)}{2X} = \sum\_{i=0}^{\frac{d + 1}{2} - 1} c\_{2i + 1} X^{2i} \\)

The crucial step in FRI is the prover's construction of a new "folded" polynomial with respect to \\(Y = X^2\\): 

\\[f^{*}(Y) = f_{E}(Y) + \alpha \cdot f_{O}(Y)\\]

Here, \\(\alpha\\) is a random value provided by the verifier. Notice that while \\(f(X)\\) has degree \\(d\\), **the new polynomial \\(f^{*}(Y)\\) has degree approximately \\(\frac{d}{2}\\)**. This degree reduction is the engine of FRI's efficiency.

Let \\(D \subset \mathbb{F}_{p} \\) be a multiplicative subgroup of even order \\(N\\), and let \\(\omega\\) be a generator of \\(D\\). The prover commits to the evaluations of \\(f(x)\\) at all elements of \\(D\\), i.e., the codewords \\(\\{f(\omega^{i})\\}\_{i=0}^{N - 1}\\).

For any \\(i \in \{0, 1, \dots \frac{N}{2}-1\}\\), we have the following relationship:

\begin{align*}
f^{\*}(w^{2i}) &= \frac{f(\omega^{i}) + f(- \omega^{i})}{2} + \alpha \cdot \frac{f(\omega^{i}) - f(- \omega^{i})}{2\omega^{i}} \\\\
               &= \frac{1}{2} \\{(1 + \alpha \omega^{-i}) \cdot f(\omega^{i}) + (1 - \alpha \omega^{-i}) \cdot f(- \omega^{i})\\} 
\end{align*}

Since \\(\omega\\) has order \\(N\\), we have that \\(\omega^{\frac{N}{2}} = -1\\) and \\(f(- \omega^{i}) = f(\omega^{\frac{N}{2} + i})\\). Thus, we can rewrite the above as

\\[f^{*}(\omega^{2i}) = \frac{1}{2} \left( (1 + \alpha \omega^{-i}) \cdot f(\omega^{i}) + (1 - \alpha \omega^{-i}) \cdot f(\omega^{\frac{N}{2} + i}) \right)\\]

In addition, consider the following three points:

- A: \\((\omega^{i}, f(\omega^{i}))\\)
- B: \\((\omega^{\frac{N}{2} + i}, f(\omega^{\frac{N}{2} + i}))\\)
- C: \\((\alpha, f^{*}(\omega^{2i}))\\)

**These three points lie on a straight line if and only if \\(f^{*}(\omega^{2i})\\) is correctly computed** from \\(f(\omega^{i})\\) and \\(f(\omega^{\frac{N}{2} + i})\\). We can see it by applying the Lagrance interpolation for A and B:

\begin{align*}
y &= f(\omega^{i}) \cdot \frac{x - \omega^{\frac{N}{2} + i}}{\omega^{i} - \omega^{\frac{N}{2} + i}} + f(\omega^{\frac{N}{2} + i}) \cdot \frac{x - \omega^{i}}{\omega^{\frac{N}{2} + i} - \omega^{i}} \\\\
  &= f(\omega^{i}) \cdot \frac{x + \omega^{i}}{\omega^{i} + \omega^{i}} + f(\omega^{\frac{N}{2} + i}) \cdot \frac{x - \omega^{i}}{-\omega^{i} - \omega^{i}} \\\\
  &= \frac{1}{2} \left( (1 + x \omega^{-i}) f(\omega^{i}) + (1 - x \omega^{-i}) f(\omega^{\frac{N}{2} + i}) \right)
\end{align*}

## FRI Protocol (Interactive)

The FRI protocol proceeds in rounds:

1. **Commitment Phase**
   - The prover computes the evaluations of the initial polynomial \\(f(X)\\) over the subgroup \\(D\\) of size \\(N\\), forming the codeword \\(\\{f(\omega^{i})\\}\_{i=0}^{N - 1}\\).
   - The verifier sends a random value \\(\alpha\\) to the prover.
   - The prover computes the evaluations of the folded polynomial \\(f^{\*}(x)\\) over the subgroup \\(D^2\\), forming the other codeword \\(\\{f^{*}(w^{2i})\\}\_{i=0}^{\frac{N}{2} - 1}\\)
   - The prover commits to these codewords by sending their Merkle roots to the verifier. The prover cannot change codewords after this commitment.
  
2. **Interactive Query Phase**
   - The verifier randomly samples an index \\(i \leftarrow \\{0, \dots, \frac{N}{2} - 1\\}\\).
   - The prover reveals the values \\(f(\omega^{i})\\), \\(f(\omega^{\frac{N}{2} + i})\\), and \\(f^{*}(\omega^{2i})\\) along with their Merkle authentication paths.
   - The verifier checks the Merkle paths to ensure the revealed values are consistent with the previously sent Merkle roots.
   - The verifier performs the collinearity check: verifies if the point \\((\alpha, f^{*}(\omega^{2i}))\\) lies on the line defined by \\((\omega^{i}, f(\omega^{i}))\\) and \\((\omega^{\frac{N}{2} + i}, f(\omega^{\frac{N}{2} + i}))\\).

3. **Degree Reduction**
   - If the collinearity check passes, the verifier is convinced that the folded polynomial \\(f^{\*}(X)\\) was correctly derived. The protocol then repeats with \\(f^{\*}(X)\\) and the subgroup \\(D^2\\). The size of the subgroup and the effective degree of the polynomial are roughly halved in each round.

This process continues for \\(\log_2{d}\\) rounds. At the final round, the remaining polynomial is a constant function, and the verifier can directly query all the remaining evaluation points to verify the degree bound. Note that with each successive round \\(j\\), the size of the domain \\(D^{j}\\) is halved. If \\(d = 2^{n}\\) and \\(N = 8d = 2^{3+n}\\), the size of the domain after \\(\log_2{d} = n\\) rounds is only \\(8\\).

## FRI Protocol (Non-Interactive)

TBD