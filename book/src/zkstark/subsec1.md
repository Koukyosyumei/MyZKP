# FRI

The goal of Fast Reed-Solomon IOP of Proximity (FRI) is to prove that a polynomial has a bounded degree.

Let \\(f(X) = \sum^{d}_{i=0} c_i X^{i}\\) be the polynomial, where \\(d\\) is the maximum degree. Following the dived-and-conquer strategy of the fast Fourie transform, this polynomial is divided into even and odd terms:

\\[f(X) = f_E(X^2) + X \cdot f_O(X^2)\\]

, where 

\begin{align*}
f_E(X^2) = \frac{f(X) + f(-X)}{2} = \sum^{\frac{d + 1}{2} - 1}_{i=0} c_{2i} X^{2i} \\\\
f_O(X^2) = \frac{f(X) - f(-X)}{2X} = \sum^{\frac{d + 1}{2} - 1}_{i=0} c_{2i + 1} X^{2i}
\end{align*}

We also introduce a new polynomial \\(f^{*}(X) = f_{E}(X) + \alpha \cdot f_{O}(X)\\), where \\(\alpha\\) is a random value.

Let \\(D \subset \mathbb{F}_{p} \ \{0\} \\) be a subgroup of even order \\(N\\) of the multiplicative group of the field, and let \\(omega\\) be the generator of \\(D\\).

Then for \\(i \in \{0, 1, \dots \frac{N}{2}-1\}\\), we have the following relationship:

\\[f^{*}(w^{2i}) = \frac{f(\omega^{i}) + f(- \omega^{i})}{2} + \alpha \cdot \frac{f(\omega^{i}) - f(- \omega^{i})}{2\omega^{i}} = \frac{1}{2} \{(1 + \alpha \omega^{-i}) \cdot f(\omega^{i}) + (1 - \alpha \omega^{-i}) \cdot f(- \omega^{i}) \}\\]

Note that since the order of \\(\omega\\) is \\(N\\), we have \\(\omega^{\frac{N}{2}} = -1\\), and \\(f(- \omega^{i}) = f(\omega^{\frac{N}{2} + i})\\). 

In addition, consider the following three points:

- A: \\((\omega^{i}, f(\omega^{i}))\\)
- B: \\((\omega^{\frac{N}{2} + i}, f(\omega^{\frac{N}{2} + i}))\\)
- C: \\((\alpha, f^{*}(\omega^{2i}))\\)

Tha magic is that those three points lie on a straight line. We can see it by applying the Lagrance interpolation for A and B:

\begin{align*}
y = f(\omega^{i}) \cdot \frac{x - \omega^{\frac{N}{2} + i}}{\omega^{i} - \omega^{\frac{N}{2} + i}} + f(\omega^{\frac{N}{2} + i}) \cdot \frac{x - \omega^{i}}{\omega^{\frac{N}{2} + i} - \omega^{i}}
\end{align*}

FRI leverage thoese relationship. Key observation from the above is as follows:

- \\(f^{*}(W^{2i})\\) is essentiaally a function of \\(f(\omega^{i})\\) and \\(f(\omega^{\frac{N}{2} + i})\\)
- While the degree of \\(f(X) = f_E{X^2} + f_O(X^2)\\) is \\(d\\), the degree of \\(f^{*}(X)\\) is \\(\frac{d}{2}\\), which is half of \\(f(X)\\).

Let's first see the overview procedures of FRI:

- The prover first computes codewords; \\(\{f(\omega^{i})\}^{N - 1}_{i=0}\\), \\(\{f_E(w^{2i})\}^{\frac{N}{2} - 1}_{i=0}\\), \\(\{f_O(w^{2i})\}^{\frac{N}{2} - 1}_{i=0}\\), and  \\(\{f^{*}(w^{2i})\}^{\frac{N}{2} - 1}_{i=0}\\).
- The prover submits Merkle roots of thoese codewords, respectively.
- The verifier randomly samples an index \\(i \rightarrow \{0, \dots, \frac{N}{2} - 1\}\\).
- The prover sends \\((f(\omega^{i}, f(\omega^{\frac{N}{2}} + i), f^{*}(\omega^{2i})))\\) along with their Merkle authentication paths.
- The verifier checks the submitted Merkle authentication paths.
- The verifier executes the colinearlity check.

After those procedures, the verifier can conclude that :

- the prover has \\(N\\) codewords corresponding to \\(f\\)
- the prover has \\(N/2 - 1\\) codewords corresponding to \\(f^{*}\\) whose degree is half of \\(f\\)

If \\(f\\)'s degree is actually equal to or less than \\(d\\), repeating the above process should stop in \\(log2 (d)\\) steps. If not, the verifier can conclude that \\(f\\)'s degree is higher than \\(d\\).
