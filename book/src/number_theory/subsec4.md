# Galois Field

We will now discuss the construction of finite fields with \\(p^n\\) elements, where \\(p\\) is a prime number and \\(n\\) is a positive integer. These fields are also known as Galois fields, denoted as \\(GF(p^n)\\). It is evident that \\(\mathbb{Z}/p\mathbb{Z}\\) is isomorphic to \\(GF(p)\\).


### Definition: Irreducible Polynomial

---

*A polynomial \\(f \in (\mathbb{Z}/p\mathbb{Z})[X]\\) of degree \\(n\\) is called irreducible over \\(\mathbb{Z}/p\mathbb{Z}\\) if it cannot be factored as a product of two polynomials of lower degree in \\((\mathbb{Z}/p\mathbb{Z})[X]\\).*

---

**Example:** In \\((\mathbb{Z}/2\mathbb{Z})[X]\\):

- \\(X^2 + X + 1\\) is irreducible
- \\(X^2 + 1 = (X + 1)^2\\) is reducible

**Implementation:**

```rust
pub trait IrreduciblePoly<F: Field>: Debug + Clone + Hash {
    fn modulus() -> &'static Polynomial<F>;
}
```

### Definition: Residue Class modulo a Polynomial

To construct \\(GF(p^n)\\), we use an irreducible polynomial of degree \\(n\\) over \\(\mathbb{Z}/p\mathbb{Z}\\).

---

*For \\(f, g \in (\mathbb{Z}/p\mathbb{Z})[X]\\), the residue class of \\(g \bmod f\\) is the set of all polynomials \\(h \in (\mathbb{Z}/p\mathbb{Z})[X]\\) such that \\(h \equiv g \pmod{f}\\). This class is denoted as:*

\\[ g + f(\mathbb{Z}/p\mathbb{Z})[X] = \\{g + hf : h \in (\mathbb{Z}/p\mathbb{Z})[X]\\} \\]

---

**Example:** In \\((\mathbb{Z}/2\mathbb{Z})[X]\\), with \\(f(X) = X^2 + X + 1\\), the residue classes \\(\bmod f\\) are:

- \\(0 + f(\mathbb{Z}/2\mathbb{Z})[X] = \\{0, X^2 + X + 1, X^2 + X, X^2 + 1, X^2, X + 1, X, 1\\}\\)
- \\(1 + f(\mathbb{Z}/2\mathbb{Z})[X] = \\{1, X^2 + X, X^2 + 1, X^2, X + 1, X, 0, X^2 + X + 1\\}\\)
- \\(X + f(\mathbb{Z}/2\mathbb{Z})[X] = \\{X, X^2 + 1, X^2, X^2 + X + 1, X + 1, 1, X^2 + X, 0\\}\\)
- \\((X + 1) + f(\mathbb{Z}/2\mathbb{Z})[X] = \\{X + 1, X^2, X^2 + X + 1, X^2 + X, 1, 0, X^2 + 1, X\\}\\)

These four residue classes form \\(GF(4)\\).

**Implementation:**

```rust
impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>>
    ExtendedFieldElement<M, P>
{
    pub fn new(poly: Polynomial<FiniteFieldElement<M>>) -> Self {
        let result = Self {
            poly: poly,
            _phantom: PhantomData,
        };
        result.reduce()
    }

    fn reduce(&self) -> Self {
        Self {
            poly: &self.poly.reduce() % P::modulus(),
            _phantom: PhantomData,
        }
    }

    pub fn degree(&self) -> isize {
        P::modulus().degree()
    }

    pub fn from_base_field(value: FiniteFieldElement<M>) -> Self {
        Self::new((Polynomial { coef: vec![value] }).reduce()).reduce()
    }
}

impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>> Field
    for ExtendedFieldElement<M, P>
{
    fn inverse(&self) -> Self {
        let mut lm = Polynomial::<FiniteFieldElement<M>>::one();
        let mut hm = Polynomial::zero();
        let mut low = self.poly.clone();
        let mut high = P::modulus().clone();

        while !low.is_zero() {
            let q = &high / &low;
            let r = &high % &low;
            let nm = hm - (&lm * &q);
            high = low;
            hm = lm;
            low = r;
            lm = nm;
        }

        Self::new(hm * high.coef[0].inverse())
    }

    fn div_ref(&self, other: &Self) -> Self {
        self.mul_ref(&other.inverse())
    }
}
```

### Theorem:

---

*If \\(f \in (\mathbb{Z}/p\mathbb{Z})[X]\\) is an irreducible polynomial of degree \\(n\\), then the residue ring \\((\mathbb{Z}/p\mathbb{Z})[X]/(f)\\) is a field with \\(p^n\\) elements, isomorphic to \\(GF(p^n)\\).*

---

**Proof** (Outline) The irreducibility of \\(f\\) ensures that \\((f)\\) is a maximal ideal in \\((\mathbb{Z}/p\mathbb{Z})[X]\\), making the quotient ring a field. The number of elements is \\(p^n\\) because there are \\(p^n\\) polynomials of degree less than \\(n\\) over \\(\mathbb{Z}/p\mathbb{Z}\\).

This construction allows us to represent elements of \\(GF(p^n)\\) as polynomials of degree less than \\(n\\) over \\(\mathbb{Z}/p\mathbb{Z}\\). Addition is performed coefficient-wise modulo \\(p\\), while multiplication is performed modulo the irreducible polynomial \\(f\\).


**Example:** To construct \\(GF(8)\\), we can use the irreducible polynomial \\(f(X) = X^3 + X + 1\\) over \\(\mathbb{Z}/2\mathbb{Z}\\). The elements of \\(GF(8)\\) are represented by:
\\[ \\{0, 1, X, X+1, X^2, X^2+1, X^2+X, X^2+X+1\\} \\]
For instance, multiplication in \\(GF(8)\\):
\\[ (X^2 + 1)(X + 1) = X^3 + X^2 + X + 1 \equiv X^2 \pmod{X^3 + X + 1} \\]

**Implementation:**

```rust
impl<M: ModulusValue + 'static, P: IrreduciblePoly<FiniteFieldElement<M>>> Ring
    for ExtendedFieldElement<M, P>
{
    fn add_ref(&self, other: &Self) -> Self {
        Self::new(&self.poly + &other.poly)
    }

    fn mul_ref(&self, other: &Self) -> Self {
        Self::new(&self.poly * &other.poly)
    }

    fn sub_ref(&self, other: &Self) -> Self {
        Self::new(&self.poly - &other.poly)
    }
}
```

### Lemma: Schwartz - Zippel Lemma

---

*Let \\(\mathbb{F}\\) be a field and \\(P: F^m \rightarrow \mathbb{F}\\) and \\(Q: \mathbb{F}^m \rightarrow \mathbb{F}\\) be two distinct multivariate polynomials of total degree at most \\(n\\). For any finite subset \\(\mathbb{S} \subseteq \mathbb{F}\\), we have:*

\begin{align}
    Pr_{u \sim \mathbb{S}^{m}}[P(u) = Q(u)] \leq \frac{n}{|\mathbb{S}|}
\end{align}

*, where \\(u\\) is drawn uniformly at random from \\(\mathbb{S}^m\\).*

---

**Proof** TBD


This lemma states that if \\(\mathbb{S}\\) is sufficiently large and \\(n\\) is relatively small, the probability that the two different polynomials return the same value for a randomly chosen input is negligibly small. In other words, if we observe \\(P(u) = Q(u)\\) for a random input \\(u\\), we can conclude with high probability that \\(P\\) and \\(Q\\) are identical polynomials.
