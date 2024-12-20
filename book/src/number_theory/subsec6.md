# Pairing

### Definition: Pairing

---

*Let \\(G_1\\) and \\(G_2\\) be cyclic groups under addition, both of prime order \\(p\\), with generators \\(P\\) and \\(Q\\) respectively:*

\begin{align}
    G_1 &= \\{0, P, 2P, ..., (p-1)P\\} \\
    G_2 &= \\{0, Q, 2Q, ..., (p-1)Q\\}
\end{align}

*Let \\(G_T\\) be a cyclic group under multiplication, also of order \\(p\\). A pairing is a map \\(e: G_1 \times G_2 \rightarrow G_T\\) that satisfies the following bilinear property:*

\begin{equation}
    e(aP, bQ) = e(P, Q)^{ab}
\end{equation} *for all \\(a, b \in \mathbb{Z}_p\\).*

---

Imagine \\(G_1\\) represents length, \\(G_2\\) represents width, and \\(G_T\\) represents area. The pairing function \\(e\\) is like calculating the area: If you double the length and triple the width, the area becomes six times larger: \\(e(2P, 3Q) = e(P, Q)^{6}\\)

### Definition: The Weil Pairing

The Weil pairing is one of the bilinear pairings for elliptic curves. We begin with its formal definition.

---

*Let \\(E\\) be an elliptic curve and \\(n\\) be a positive integer. For points \\(P, Q \in E[n]\\), where \\(E[n]\\) denotes the \\(n\\)-torsion subgroup of \\(E\\), we define the Weil pairing \\(e_n(P, Q)\\) as follows:*

Let \\(f_P\\) and \\(f_Q\\) be rational functions on \\(E\\) satisfying:
\begin{align}
    div(f_P) &= n[P] - n[\mathcal{O}] \\
    div(f_Q) &= n[Q] - n[\mathcal{O}]
\end{align}

*Then, for an arbitrary point \\(S \in E\\) such that \\(S \notin \\{\mathcal{O}, P, -Q, P-Q\\}\\), the Weil pairing is given by:*

\begin{equation}
e_n(P, Q) = \frac{f_P(Q + S)}{f_P(S)} /\ \frac{f_Q(P - S)}{f_Q(-S)}
\end{equation}

---

We introduce a crucial theorem about a specific rational function on elliptic curves to construct the functions required for the Weil pairing.

### Theorem

---

*Let \\(E\\) be an elliptic curve over a field \\(K\\), and let \\(P = (x_P, y_P)\\) and \\(Q = (x_Q, y_Q)\\) be non-zero points on \\(E\\). Define \\(\lambda\\) as:*

\begin{equation}
    \lambda = \begin{cases}
        \hbox{slope of the line through \\(P\\) and \\(Q\\)} &\quad \hbox{if \\(P \neq Q\\)} \\\\
        \hbox{slope of the tangent line to \\(E\\) at \\(P\\)} &\quad \hbox{if \\(P = Q\\)} \\\\
        \infty &\quad \hbox{if the line is vertical}
    \end{cases}
\end{equation}

*Then, the function \\(g_{P,Q}: E \to K\\) defined by:*
\begin{equation}
g_{P,Q} = \begin{cases}
\frac{y - y_P - \lambda(x - x_P)}{x + x_P + x_Q - \lambda^2} &\quad \hbox{if } \lambda \neq \infty \\\\
x - x_P &\quad \hbox{if } \lambda = \infty
\end{cases}
\end{equation} *has the following divisor:*

\begin{equation}
div(g_{P,Q}) = [P] + [Q] - [P + Q] - [\mathcal{O}]
\end{equation}

---

**Proof:** We consider two cases based on the value of \\(\lambda\\).

Case 1: \\(\lambda \neq \infty\\)

Let \\(y = \lambda x + v\\) be the line through \\(P\\) and \\(Q\\) (or the tangent line at \\(P\\) if \\(P = Q\\)). This line intersects \\(E\\) at three points: \\(P\\), \\(Q\\), and \\(-P-Q\\). Thus,
\begin{equation}
div(y - \lambda x - v) = [P] + [Q] + [-P - Q] - 3[\mathcal{O}]
\end{equation}
Vertical lines intersect \\(E\\) at points and their negatives, so:
\begin{equation}
div(x - x_{P+Q}) = [P + Q] + [-P - Q] - 2[\mathcal{O}]
\end{equation}
It follows that \\(g_{P,Q} = \frac{y - \lambda x - v}{x - x_{P+Q}}\\) has the desired divisor.

Case 2: \\(\lambda = \infty\\)

In this case, \\(P + Q = \mathcal{O}\\), so we want \\(g_{P,Q}\\) to have divisor \\([P] + [-P] - 2[\mathcal{O}]\\). The function \\(x - x_P\\) has this divisor.
\end{proof}

### Theorem: Miller's Algorithm

---

*Let \\(m \geq 1\\) and write its binary expansion as:*
\begin{equation}
m = m_0 + m_1 \cdot 2 + m_2 \cdot 2^2 + \cdots + m_{n-1} \cdot 2^{n-1}
\end{equation}
*where \\(m_i \in \\{0, 1\\}\\) and \\(m_{n-1} \neq 0\\).*

*The following algorithm, using the function \\(g_{P,Q}\\) defined in the previous theorem, returns a function \\(f_P\\) whose divisor satisfies:*

\begin{equation}
    div(f_P) = m[P] - m[P] - (m - 1)[\mathcal{O}]
\end{equation}

===========================

**Miller's Algorithm**

1. Set \\(T = P\\) and \\(f = 1\\)
2. **For** \\(i \gets n-2 \cdots 0\\) **do**
3. &nbsp;&nbsp;&nbsp;&nbsp;Set \\(f = f^2 \cdot g_{T, T}\\)
4. &nbsp;&nbsp;&nbsp;&nbsp;Set \\(T = 2T\\)
5. &nbsp;&nbsp;&nbsp;&nbsp;**If** \\(m_i = 1\\) **then**
6. &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Set \\(f = f \cdot g_{T, P}\\)
7. &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Set \\(T = T + P\\)
8. &nbsp;&nbsp;&nbsp;&nbsp;**End if**
9. **End for**
10. **Return** \\(f\\)
    
===========================

---

**Proof:** TBD

**Implementation:**

```rust
pub fn get_lambda<F: Field, E: EllipticCurve>(
    p: &EllipticCurvePoint<F, E>,
    q: &EllipticCurvePoint<F, E>,
    r: &EllipticCurvePoint<F, E>,
) -> F {
    let p_x = p.x.as_ref().unwrap();
    let p_y = p.y.as_ref().unwrap();
    let q_x = q.x.as_ref().unwrap();
    // let q_y = q.y.clone().unwrap();
    let r_x = r.x.as_ref().unwrap();
    let r_y = r.y.as_ref().unwrap();

    if (p == q && *p_y == F::zero()) || (p != q && *p_x == *q_x) {
        return r_x.sub_ref(&p_x);
    }
    let slope = p.line_slope(&q);
    let numerator = (r_y.sub_ref(&p_y)).sub_ref(&slope.mul_ref(&(r_x.sub_ref(&p_x))));
    let denominator = r_x
        .add_ref(&p_x)
        .add_ref(&q_x)
        .sub_ref(&slope.mul_ref(&slope));
    return numerator / denominator;
}

pub fn miller<F: Field, E: EllipticCurve>(
    p: &EllipticCurvePoint<F, E>,
    q: &EllipticCurvePoint<F, E>,
    m: &BigInt,
) -> (F, EllipticCurvePoint<F, E>) {
    if p == q {
        return (F::one(), p.clone());
    }

    let mut f = F::one();
    let mut t = p.clone();

    for i in (0..(m.bits() - 1)).rev() {
        f = (f.mul_ref(&f)) * (get_lambda(&t, &t, &q));
        t = t.add_ref(&t);
        if m.bit(i) {
            f = f * (get_lambda(&t, &p, &q));
            t = t.add_ref(&p);
        }
    }

    (f, t)
}
```
