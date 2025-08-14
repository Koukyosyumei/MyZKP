# KZG

We consider the polynomial commitment shceme, where the target polynomial \\(f \in F_{q}[X]\\) such that \\(f(x) = \sum_{i=0}^{d} c_{i} x^{i}\\).

## Setup

Let \\(D\\) be a maximu degree such that \\(d \leq D\\).

\begin{align*}
    &\mathrm{SK} : \alpha \xleftarrow{\text{random}} F_{q} \\\\
    &\mathrm{PK} : (\\{g_1^{(\alpha^{i})}\\}^{D}_{i=0}, g_2, g^{\alpha}_2)
\end{align*}

## Commitment

\\[
    c = \Pi_{i=0}^{d} (g_1^{(\alpha^{i})})^{c_i}
\\]

We can observe that

\begin{align*}
    c &= \Pi_{i=0}^{d} (g_1^{(\alpha^{i})})^{c_i} = g_1^{\sum_{i=0}^{d} c_i \alpha^{i}} = g_1^{f(\alpha)}
\end{align*}

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