## Semigroup, Group, Ring

---

***Definition: Semigroup***

A pair \\((H, \circ)\\), where \\(H\\) is a non-empty set and \\(\circ\\) is an associative binary operation on \\(H\\), is called a semigroup.

**Example**: The set of positive integers under multiplication modulo \\(n\\) forms a semigroup. For instance, with \\(n = 6\\), the elements \\(\{1, 2, 3, 4, 5\}\\) under multiplication modulo 6 form a semigroup, since multiplication modulo 6 is associative.

---

***Definition: Abelian Semigroup***

A semigroup whose operation is commutative is called an abelian semigroup.

**Example**: The set of natural numbers under addition modulo \\(n\\) forms an abelian semigroup. For \\(n = 7\\), addition modulo 7 is both associative and commutative, so it is an abelian semigroup.

---

***Definition: Identity Element***

An element \\(e \in H\\) is an identity element of \\(H\\) if it satisfies \\(e \circ a = a \circ e = a\\) for any \\(a \in H\\).

**Example**: 0 is the identity element for addition modulo \\(n\\). For example, \\(0 + a \bmod 5 = a + 0 \bmod 5 = a\\). Similarly, 1 is the identity element for multiplication modulo \\(n\\). For example, \\(1 \times a \bmod 7 = a \times 1 \bmod 7 = a\\).

---

***Definition: Monoid***

A semigroup with an identity element is called a monoid.

**Example**: The set of non-negative integers under addition modulo \\(n\\) forms a monoid. For \\(n = 5\\), the set \\(\{0, 1, 2, 3, 4\}\\) under addition modulo 5 forms a monoid with 0 as the identity element.

---

***Definition: Inverse***

For an element \\(a \in H\\), an element \\(b \in H\\) is an inverse of \\(a\\) if \\(a \circ b = b \circ a = e\\), where \\(e\\) is the identity element.

**Example**: In modulo \\(n\\) arithmetic (addition), the inverse of an element exists if it can cancel itself out to yield the identity element. In the set of integers modulo 7, the inverse of 3 is 5, because \\(3 \times 5 \bmod 7 = 1\\), where 1 is the identity element for multiplication.

---

***Definition: Group***

A monoid in which every element has an inverse is called a group.

**Example**: The set of integers modulo a prime \\(p\\) under multiplication forms a group (Can you prove it?). For instance, in \\(\mathbb{Z}/5\mathbb{Z}\\), every non-zero element \\(\{1 + 5\mathbb{Z}, 2 + 5\mathbb{Z}, 3 + 5\mathbb{Z}, 4 + 5\mathbb{Z}\}\\) has an inverse, making it a group.

---

***Definition: Order of a Group***

The order of a group is the number of elements in the group.

**Example**: The group of integers modulo 4 under addition has order 4, because the set of elements is \\(\{0, 1, 2, 3\}\\).

---

***Definition: Ring***

A triple \\((R, +, \cdot)\\) is a ring if \\((R, +)\\) is an abelian group, \\((R, \cdot)\\) is a semigroup, and the distributive property holds: \\(x \cdot (y + z) = (x \cdot y) + (x \cdot z)\\) and \\((x + y) \cdot z = (x \cdot z) + (y \cdot z)\\) for all \\(x, y, z \in R\\).

**Example**: The set of integers with usual addition and multiplication modulo \\(n\\) forms a ring. For example, in \\(\mathbb{Z}/6\mathbb{Z}\\), addition and multiplication modulo 6 form a ring.

```rust
use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

pub trait Ring:
    Sized
    + Clone
    + PartialEq
    + fmt::Display
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Mul<i64, Output = Self>
    + Neg<Output = Self>
{
    // A ring is an algebraic structure with addition and multiplication
    fn zero() -> Self;
    fn one() -> Self;
}
```

---

***Definition: Communicative Ring***

A ring is called a commutative ring if its multiplication operation is commutative.

**Example** The set of real numbers under usual addition and multiplication forms a commutative ring.

---

***Definition: Field***

A commutative ring with a multiplicative identity element where every non-zero element has a multiplicative inverse is called a field.

**Example** The set of rational numbers under usual addition and multiplication forms a field.

```rust
use crate::modules::ring::Ring;

pub trait Field: Ring + Div<Output = Self> {
    fn inverse(&self) -> Self;
}
```

---


