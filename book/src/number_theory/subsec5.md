# Elliptic curve

### Definition: Elliptic Curve
    
---

*An elliptic curve \\(E\\) over a finite field \\(\mathbb{F}\_{p}\\) is defined by the equation:*

\begin{equation}
    y^2 = x^3 + a x + b
\end{equation}
    
*, where \\(a, b \in \mathbb{F}\_{p}\\) and the discriminant \\(\Delta_{E} = 4a^3 + 27b^2 \neq 0\\).*

---

```rust
use crate::modules::field::Field;

#[derive(Debug, Clone, PartialEq)]
pub struct EllipticCurve<F: Field> {
    pub a: F,
    pub b: F,
}
```

### Definition: \\(\mathbb{F}\_{p}\\)-Rational Point

---

*An \\(\mathbb{F}\_{p}\\)-rational point on an elliptic curve \\(E\\) is a point \\((x, y)\\) where both \\(x\\) and \\(y\\) are elements of \\(\mathbb{F}\_{p}\\) and satisfy the curve equation, or the point at infinity \\(\mathcal{O}\\).*

---

Example: For \\(E: y^2 = x^3 + 3x + 4\\) over \\(\mathbb{F}\_7\\), the \\(\mathbb{F}\_{7}\\)-rational points are:
\\(\\{(0, 2), (0, 5), (1, 1), (1, 6), (2, 2), (2, 5), (5, 1), (5, 6), \mathcal{O}\\}\\)

### Definition: The point at infinity

---

The point at infinity denoted \\(\mathcal{O}\\), is a special point on the elliptic curve that serves as the identity element for the group operation. It can be visualized as the point where all vertical lines on the curve meet.

---

```rust
#[derive(Debug, Clone, PartialEq)]
pub struct EllipticCurvePoint<F: Field> {
    pub x: Option<F>,
    pub y: Option<F>,
    pub curve: EllipticCurve<F>,
}
```

### Definition: Addition on Elliptic Curve

---

*For an elliptic curve \\(E: y^2 = x^3 + ax + b\\), the addition of points \\(P\\) and \\(Q\\) to get \\(R = P + Q\\) is defined as follows:*

   
   - *If \\(P = \mathcal{O}\\), \\(R = Q\\).*
   - *If \\(Q = \mathcal{O}\\), \\(R = P\\).*
   - *Otherwise, let \\(P = (x_P, y_P), Q = (x_Q, y_Q)\\). Then:*
       
       - *If \\(y_P = -y_Q\\), \\(R = \mathcal{O}\\)*
       - *If \\(y_P \neq -y_Q\\), \\(R = (x_R = \lambda^2 - x_P - x_Q, y_R = \lambda(x_P - x_R) - y_P)\\), where \\(\lambda = \\)* \\( \begin{cases}
               \frac{y_P - y_Q}{x_P - x_Q} \quad &\hbox{If } (x_P \neq x_Q) \\\\
               \frac{3^{2}_{P} + a}{2y_P} \quad &\hbox{Otherwise}
           \end{cases}\\)
       
---


Example: On \\(E: y^2 = x^3 + 2x + 3\\) over \\(\mathbb{F}_{7}\\), let \\(P = (5, 1)\\) and \\(Q = (4, 4)\\). Then, \\(P + Q = (0, 5)\\), where \\(\lambda = \frac{1 - 4}{5 - 4} \equiv 4 \bmod 7\\).

```rust
impl<F: Field> EllipticCurvePoint<F> {
    fn new(x: F, y: F, curve: EllipticCurve<F>) -> Self {
        EllipticCurvePoint {
            x: Some(x),
            y: Some(y),
            curve: curve,
        }
    }

    pub fn point_at_infinity(curve: EllipticCurve<F>) -> Self {
        EllipticCurvePoint {
            x: None,
            y: None,
            curve: curve,
        }
    }

    pub fn is_point_at_infinity(&self) -> bool {
        self.x.is_none() || self.y.is_none()
    }

    pub fn line_slope(&self, other: Self) -> F {
        let x1 = self.x.clone().unwrap();
        let y1 = self.y.clone().unwrap();
        let x2 = other.x.clone().unwrap();
        let y2 = other.y.clone().unwrap();

        if self.x.clone() == other.x.clone() {
            ((x1.clone() * x1.clone()) * (3_i64) + self.curve.a.clone()) / (y1.clone() * (2_i64))
        } else {
            (y2.clone() - y1.clone()) / (x2.clone() - x1.clone())
        }
    }
}

impl<F: Field> Add for EllipticCurvePoint<F> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        if self.is_point_at_infinity() {
            return other;
        }
        if other.is_point_at_infinity() {
            return self;
        }

        let m = self.line_slope(other.clone());

        if self.x == other.x {
            if self.y != other.y {
                return EllipticCurvePoint::point_at_infinity(self.curve.clone());
            } else {
                let x1 = self.x.clone().unwrap();
                let y1 = self.y.clone().unwrap();

                let x3 = m.clone() * m.clone() - x1.clone() - x1.clone();
                let y3 = m * (x1 - x3.clone()) - y1;

                return EllipticCurvePoint::new(x3, y3, self.curve.clone());
            }
        } else {
            let x1 = self.x.clone().unwrap();
            let y1 = self.y.clone().unwrap();
            let x2 = other.x.clone().unwrap();
            let x3 = m.clone() * m.clone() - x1.clone() - x2.clone();
            let y3 = m * (x1 - x3.clone()) - y1;

            return EllipticCurvePoint::new(x3, y3, self.curve.clone());
        }
    }
}
```

### Definition: Mordell-Weil Group

---

*The Mordell-Weil group of an elliptic curve \\(E\\) is the group of rational points on \\(E\\) under the addition operation defined above.*

---

Example: For \\(E: y^2 = x^3 + x + 6\\) over \\(\mathbb{F}\_{11}\\), the Mordell-Weil group is the set of all \\(\mathbb{F}\_{11}\\)-rational points on \\(E\\) with the elliptic curve addition operation.

### Definition: Group Order

---

*The group order of an elliptic curve \\(E\\) over \\(\mathbb{F}\_{p}\\), denoted \\(\\#E(\mathbb{F}\_{p})\\), is the number of \\(\mathbb{F}\_{p}\\)-rational points on \\(E\\), including the point at infinity.*

---

Example: For \\(E: y^2 = x^3 + x + 6\\) over \\(\mathbb{F}\_{11}\\), \\(\\#E(\mathbb{F}\_{11}) = 13\\).

### Theorem: Hasse-Weil

---

*Let \\(\\#E(\mathbb{F}\_{p})\\) be the group order of the elliptic curve \\(E\\) over \\(\mathbb{F}\_{p}\\). Then:*

\begin{equation}
    p + 1 - 2 \sqrt{p} \leq \\#E \leq p + 1 + 2 \sqrt{p}
\end{equation}

---

**Example:** For an elliptic curve over \\(\mathbb{F}_{23}\\), the Hasse-Weil theorem guarantees that: 

\begin{equation*}
23 + 1 - 2 \sqrt{23} \simeq 14.42 \\#E(\mathbb{F}_{23}) \geq 23 + 1 + 2 \sqrt{23} \simeq 33.58
\end{equation*}

### Definition: Point Order


The order of a point \\(P\\) on an elliptic curve is the smallest positive integer \\(n\\) such that \\(nP = \mathcal{O}\\) (where \\(nP\\) denotes \\(P\\) added to itself \\(n\\) times). We also denote the set of points of order \\(n\\), also called \textit{torsion} group, by 

\begin{equation}
    E[n] = \\{P \in E: [n]P = \mathcal{O}\\}
\end{equation}


Example: On \\(E: y^2 = x^3 + 2x + 2\\) over \\(\mathbb{F}_{17}\\), the point \\(P = (5, 1)\\) has order 18 because \\(18P = \mathcal{O}\\), and no smaller positive multiple of \\(P\\) equals \\(\mathcal{O}\\).

The intuitive view is that if you continue to add points to themselves (doubling, tripling, etc.), the lines drawn between the prior point and the next point will eventually become more vertical. When the line becomes vertical, it does not intersect the elliptic curve at any finite point. In elliptic curve geometry, a vertical line is considered to "intersect" the curve at a special place called the "point at infinity," This point is like a north pole in geographic terms: no matter which direction you go, if the line becomes vertical (reaching infinitely high), it converges to this point.

### Definition: Field Extension

---

*Let \\(F\\) and \\(L\\) be fields. If \\(F \subseteq L\\) and the operations of \\(F\\) are the same as those of \\(L\\), we call \\(L\\) a field extension of \\(F\\). This is denoted as \\(L/F\\).*

---

A field extension \\(L/F\\) naturally gives \\(L\\) the structure of a vector space over \\(F\\). The dimension of this vector space is called the degree of the extension.

**Examples:**

- \\(\mathbb{C}/\mathbb{R}\\) is a field extension with basis \\(\\{1, i\\}\\) and degree 2.
- \\(\mathbb{R}/\mathbb{Q}\\) is an infinite degree extension.
- \\(\mathbb{Q}(\sqrt{2})/\mathbb{Q}\\) is a degree 2 extension with basis \\(\\{1, \sqrt{2}\\}\\).


### Definition: Algebraic Extension

---

*A field extension \\(L/K\\) is called algebraic if every element \\(\alpha \in L\\) is algebraic over \\(K\\), i.e., \\(\alpha\\) is the root of some non-zero polynomial with coefficients in \\(K\\).*

*For an algebraic element \\(\alpha \in L\\) over \\(K\\), we denote by \\(K(\alpha)\\) the smallest field containing both \\(K\\) and \\(\alpha\\).*

---

**Example:**

- \\(\mathbb{C}/\mathbb{R}\\) is algebraic: any \\(z = a + bi \in \mathbb{C}\\) is a root of \\(x^2 - 2ax + (a^2 + b^2) \in \mathbb{R}[x]\\).
- \\(\mathbb{Q}(\sqrt[3]{2})/\mathbb{Q}\\) is algebraic: \\(\sqrt[3]{2}\\) is a root of \\(x^3 - 2 \in \mathbb{Q}[x]\\).
- \\(\mathbb{R}/\mathbb{Q}\\) is not algebraic (a field extension that is not algebraic is called \textit{transcendental}).


### Definition: Field of Rational Functions

---

*Let \\(K\\) be a field and \\(X\\) be indeterminate. The field of rational functions over \\(K\\), denoted \\(K(X)\\), is defined as:*

\begin{equation*}
    K(X) = \left\\{ \frac{f(X)}{g(X)} \middle|\ f(X), g(X) \in K[X], g(X) \neq 0 \right\\}
\end{equation*}

*where \\(K[X]\\) is the ring of polynomials in \\(X\\) with coefficients in \\(K\\).*

---

\\(K(X)\\) can be viewed as the field obtained by adjoining a letter \\(X\\) to \\(K\\). This construction generalizes to multiple variables, e.g., \\(K(X,Y)\\).

The concept of a function field naturally arises in the context of algebraic curves, particularly elliptic curves. Intuitively, the function field encapsulates the algebraic structure of rational functions on the curve.

We first construct the coordinate ring for an elliptic curve \\(E: y^2 = x^3 + ax + b\\). Consider functions \\(X: E \to K\\) and \\(Y: E \to K\\) that extract the \\(x\\) and \\(y\\) coordinates, respectively, from an arbitrary point \\(P \in E\\). These functions generate the polynomial ring \\(K[X,Y]\\), subject to the relation \\(Y^2 = X^3 + aX + b\\).

To put it simply, the function field is a field that consists of all functions that determine the value based on the point on the curve.

### Definition: Coordinate Ring of an Elliptic Curve

---

*The **coordinate ring** of an elliptic curve \\(E: y^2 = x^3 + ax + b\\) over a field \\(K\\) is defined as:*
\begin{equation}
    K[E] = K[X, Y]/(Y^2 - X^3 - aX - b)
\end{equation}

---

In other words, we can view \\(K[E]\\) as a ring representing all polynomial functions on \\(E\\). Recall that \\(K[X, Y]\\) is the polynomial ring in two variables \\(X\\) and \\(Y\\) over the field K, meaning that it contains all polynomials in \\(X\\) and \\(Y\\) with coefficients from \\(K\\). Then, the notation \\(K[X, Y]/(Y^2 - X^3 - aX - b)\\) denotes the quotient ring obtained by taking \\(K[X, Y]\\) and "modding out" by the ideal \\((Y^2 - X^3 - aX - b)\\).

**Example** For example, for an elliptic curve \\(E: y^2 = x^3 - x\\) over \\(\mathbb{Q}\\), some elements of the coordinate ring \\(\mathbb{Q}[E]\\) include:


- Constants: \\(3, -2, \frac{1}{7}, \ldots\\)
- Linear functions: \\(X, Y, 2X+3Y, \ldots\\)
- Quadratic functions: \\(X^2, XY, Y^2 (= X^3 - X), \ldots\\)
- Higher-degree functions: \\(X^3, X^2Y, XY^2 (= X^4 - X^2), \ldots\\)

### Definition: Function Field of an Elliptic Curve

Then, the function field is defined as follows:

---

*Let \\(E: y^2 = x^3 + ax + b\\) be an elliptic curve over a field \\(K\\). The **function field** of \\(E\\), denoted \\(K(E)\\), is defined as:*
\begin{equation}
    K(E) = \left\\{\frac{f}{g} \,\middle|\, f, g \in K[E], g \neq 0 \right\\}
\end{equation}
*where \\(K[E] = K[X, Y]/(Y^2 - X^3 - aX - b)\\) is the coordinate ring of \\(E\\).*

---

\\(K(E)\\) can be viewed as the field of all rational functions on \\(E\\).

### Definition: Zero/Pole of a Function

---

*Let \\(h \in K(E)\\) be a non-zero function. A point \\(P \in E\\) is called a **zero** of \\(h\\) if \\(h(P) = 0\\).*

---

---

*Let \\(h \in K(E)\\) be a non-zero function. A point \\(P \in E\\) is called a **pole** of \\(h\\) if \\(h\\) is not defined at \\(P\\) or, equivalently if \\(1/h\\) has a zero at \\(P\\).*

---

**Example** Consider the elliptic curve \\(E: Y^2 = X^3 - X\\) over a field \\(K\\) of characteristic \\(\neq 2, 3\\). Let \\(P_{-1} = (-1, 0)\\), \\(P_0 = (0, 0)\\), and \\(P_1 = (1, 0)\\) be the points where \\(Y = 0\\).


- The function \\(Y \in K(E)\\) has three simple zeros: \\(P_{-1}\\), \\(P_0\\), and \\(P_1\\).
- The function \\(X \in K(E)\\) has a double zero at \\(P_0\\) (since \\(P_0 = -P_0\\)).
- The function \\(X - 1 \in K(E)\\) has a simple zero at \\(P_1\\).
- The function \\(X^2 - 1 \in K(E)\\) has two simple zeros at \\(P_{-1}\\) and \\(P_1\\).
- The function \\(\frac{Y}{X} \in K(E)\\) has a simple zero at \\(P_{-1}\\) and a simple pole at \\(P_0\\).


An important property of functions in \\(K(E)\\) is that they have the same number of zeros and poles when counted with multiplicity. This is a consequence of a more general result known as the Degree-Genus Formula.

### Theorem: Degree-Genus Formula for Elliptic Curves

---

*Let \\(f \in K(E)\\) be a non-zero rational function on an elliptic curve \\(E\\). Then:*

\begin{equation}
    \sum_{P \in E} ord_{P(f)} = 0
\end{equation}

*, where \\(ord_{P(f)}\\) denotes the order of \\(f\\) at \\(P\\), which is positive for zeros and negative for poles.*

---

This theorem implies that the total number of zeros (counting multiplicity) equals the total number of poles for any non-zero function in \\(K(E)\\).

We now introduce a powerful tool for analyzing functions on elliptic curves: the concept of divisors.

### Definition: Divisor of a Function on an Elliptic Curve

---

*Let \\(E: Y^2 = X^3 + AX + B\\) be an elliptic curve over a field \\(K\\), and let \\(f \in K(E)\\) be a non-zero rational function on \\(E\\). The **divisor** of \\(f\\), denoted \\(div(f)\\), is defined as:*
\begin{equation}
    div(f) = \sum_{P \in E} ord_P(f) [P]
\end{equation}
*, where \\(ord_P(f)\\) is the order of \\(f\\) at \\(P\\) (positive for zeros, negative for poles), and the sum is taken over all points \\(P \in E\\), including the point at infinity. This sum has only finitely many non-zero terms.*


---

Note that this \\(\sum_{P \in E}\\) is a symbolic summation, and we do not calculate the concrete numerical value of a divisor.

**Example** Consider the elliptic curve \\(E: Y^2 = X^3 - X\\) over \\(\mathbb{Q}\\).

- For \\(f = X\\), we have \\(div(X) = 2[(0,0)] - 2[\mathcal{O}]\\).
- For \\(g = Y\\), we have \\(div(Y) = [(1,0)] + [(0,0)] + [(-1,0)] - 3[\mathcal{O}]\\).
- For \\(h = \frac{X-1}{Y}\\), we have \\(div(h) = [(1,0)] - [(-1,0)]\\).


### Definition: Divisor on an Elliptic Curve

The concept of divisors can be extended to the elliptic curve itself:

---

*A **divisor** \\(D\\) on an elliptic curve \\(E\\) is a formal sum*
\begin{equation}
    D = \sum_{P \in E} n_P [P]
\end{equation}
*where \\(n_P \in \mathbb{Z}\\) and \\(n_P = 0\\) for all but finitely many \\(P\\).*

---

**Example** On the curve \\(E: Y^2 = X^3 - X\\):

- \\(D_1 = 3[(0,0)] - 2[(1,1)] - [(2,\sqrt{6})]\\) is a divisor.
- \\(D_2 = [(1,0)] + [(-1,0)] - 2[\mathcal{O}]\\) is a divisor.
- \\(D_3 = \sum_{P \in E[2]} [P] - 4[\mathcal{O}]\\) is a divisor (where \\(E[2]\\) are the 2-torsion points).


To quantify the properties of divisors, we introduce two important metrics:


### Definition: Degree/Sum of a Divisor

---

*The **degree** of a divisor \\(D = \sum_{P \in E} n_P [P]\\) is defined as:*
\begin{equation}
    deg(D) = \sum_{P \in E} n_P
\end{equation}

---

---

*The **sum** of a divisor \\(D = \sum_{P \in E} n_P [P]\\) is defined as:*
\begin{equation}
    Sum(D) = \sum_{P \in E} n_P P
\end{equation}
*where \\(n_P P\\) denotes the point addition of \\(P\\) to itself \\(n_P\\) times in the group law of \\(E\\).*

---

**Example** For the divisors in the previous example:

- \\(deg(D_1) = 3 - 2 - 1 = 0\\)
- \\(deg(D_2) = 1 + 1 - 2 = 0\\)
- \\(deg(D_3) = 4 - 4 = 0\\)
- \\(Sum(D_2) = (1,0) + (-1,0) - 2\mathcal{O} = \mathcal{O}\\) (since \\((1,0)\\) and \\((-1,0)\\) are 2-torsion points)

### Theorem 

The following theorem characterizes divisors of functions and provides a criterion for when a divisor is the divisor of a function:

---

*Let \\(E\\) be an elliptic curve over a field \\(K\\).*


- *If \\(f, f' \in K(E)\\) are non-zero rational functions on \\(E\\) with \\(div(f) = div(f')\\), then there exists a non-zero constant \\(c \in K^\*\\) such that \\(f = cf'\\).*
- *A divisor \\(D\\) on \\(E\\) is the divisor of a rational function on \\(E\\) if and only if \\(deg(D) = 0\\) and \\(Sum(D) = \mathcal{O}\\).*


---

**Example** On \\(E: Y^2 = X^3 - X\\):

- The function \\(f = \frac{Y}{X}\\) has \\(div(f) = [(1,0)] + [(-1,0)] - 2[(0,0)]\\). Note that \\(deg(div(f)) = 0\\) and \\(Sum(div(f)) = (1,0) + (-1,0) - 2(0,0) = \mathcal{O}\\).
- The divisor \\(D = 2[(1,1)] - [(2,\sqrt{6})] - [(0,0)]\\) has \\(deg(D) = 0\\), but \\(Sum(D) \neq \mathcal{O}\\). Therefore, \\(D\\) is not the divisor of any rational function on \\(E\\).