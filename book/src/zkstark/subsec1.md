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