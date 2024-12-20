# Semigroup, Group, Ring

### Definition: Semigroup

---

*A pair \\((H, \circ)\\), where \\(H\\) is a non-empty set and \\(\circ\\) is an associative binary operation on \\(H\\), is called a semigroup.*

---

**Example:** The set of positive integers under multiplication modulo \\(n\\) forms a semigroup. For instance, with \\(n = 6\\), the elements \\(\\{1, 2, 3, 4, 5\\}\\) under multiplication modulo 6 form a semigroup, since multiplication modulo 6 is associative.


### Definition: Abelian Semigroup

---

*A semigroup whose operation is commutative is called an abelian semigroup.*

---

**Example:** The set of natural numbers under addition modulo \\(n\\) forms an abelian semigroup. For \\(n = 7\\), addition modulo 7 is both associative and commutative, so it is an abelian semigroup.


### Definition: Identity Element

---

*An element \\(e \in H\\) is an identity element of \\(H\\) if it satisfies \\(e \circ a = a \circ e = a\\) for any \\(a \in H\\).*

---

**Example:** 0 is the identity element for addition modulo \\(n\\). For example, \\(0 + a \bmod 5 = a + 0 \bmod 5 = a\\). Similarly, 1 is the identity element for multiplication modulo \\(n\\). For example, \\(1 \times a \bmod 7 = a \times 1 \bmod 7 = a\\).

### Definition: Monoid

---

*A semigroup with an identity element is called a monoid.*

---

**Example:** The set of non-negative integers under addition modulo \\(n\\) forms a monoid. For \\(n = 5\\), the set \\(\\{0, 1, 2, 3, 4\\}\\) under addition modulo 5 forms a monoid with 0 as the identity element.

### Definition: Inverse

---

*For an element \\(a \in H\\), an element \\(b \in H\\) is an inverse of \\(a\\) if \\(a \circ b = b \circ a = e\\), where \\(e\\) is the identity element.*

---

**Example:** In modulo \\(n\\) arithmetic (addition), the inverse of an element exists if it can cancel itself out to yield the identity element. In the set of integers modulo 7, the inverse of 3 is 5, because \\(3 \times 5 \bmod 7 = 1\\), where 1 is the identity element for multiplication.

### Definition: Group

---

*A monoid in which every element has an inverse is called a group.*

---

**Example:** The set of integers modulo a prime \\(p\\) under multiplication forms a group (Can you prove it?). For instance, in \\(\mathbb{Z}/5\mathbb{Z}\\), every non-zero element \\(\\{1 + 5\mathbb{Z}, 2 + 5\mathbb{Z}, 3 + 5\mathbb{Z}, 4 + 5\mathbb{Z}\\}\\) has an inverse, making it a group.

### Definition: Order of a Group

---

*The order of a group is the number of elements in the group.*

---

**Example:** The group of integers modulo 4 under addition has order 4, because the set of elements is \\(\\{0, 1, 2, 3\\}\\).


### Definition: Ring

---

*A triple \\((R, +, \cdot)\\) is a ring if \\((R, +)\\) is an abelian group, \\((R, \cdot)\\) is a semigroup, and the distributive property holds: \\(x \cdot (y + z) = (x \cdot y) + (x \cdot z)\\) and \\((x + y) \cdot z = (x \cdot z) + (y \cdot z)\\) for all \\(x, y, z \in R\\).*

---

**Example:** The set of integers with usual addition and multiplication modulo \\(n\\) forms a ring. For example, in \\(\mathbb{Z}/6\mathbb{Z}\\), addition and multiplication modulo 6 form a ring.

**Implementation:**

```rust
use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

pub trait Ring:
    Sized
    + Clone
    + PartialEq
    + fmt::Display
    + Add<Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + Sub<Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + Mul<Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + Neg<Output = Self>
    + One
    + Zero
{
    // A ring is an algebraic structure with addition and multiplication
    fn add_ref(&self, rhs: &Self) -> Self;
    fn sub_ref(&self, rhs: &Self) -> Self;
    fn mul_ref(&self, rhs: &Self) -> Self;

    // Utility functions
    fn pow<M: Into<BigInt>>(&self, n: M) -> Self;
    fn get_value(&self) -> BigInt;
    fn from_value<M: Into<BigInt>>(value: M) -> Self;
    fn random_element(exclude_elements: &[Self]) -> Self;
}
```

### Definition: Communicative Ring

---

*A ring is called a commutative ring if its multiplication operation is commutative.*

---

**Example** The set of real numbers under usual addition and multiplication forms a commutative ring.

### Definition: Field

---

*A commutative ring with a multiplicative identity element where every non-zero element has a multiplicative inverse is called a field.*

---

**Example** The set of rational numbers under usual addition and multiplication forms a field.

**Implementation:**

```rust
pub trait Field: Ring + Div<Output = Self> + PartialEq + Eq + Hash {
    /// Computes the multiplicative inverse of the element.
    fn inverse(&self) -> Self;
    fn div_ref(&self, other: &Self) -> Self;
}
```


### Definition: Residue Class

---

*The residue class of \\( a \\) modulo \\( m \\), denoted as \\( a + m\mathbb{Z} \\), is the set \\( \\{b : b \equiv a \pmod{m}\\\} \\).*

---

**Example:** For \\( m = 3 \\), the residue class of 2 is \\( 2 + 3\mathbb{Z} = \\{\ldots, -4, -1, 2, 5, 8, \ldots\\} \\).

### Definition: Inverse of Residue Class

---

*We denote the set of all residue classes modulo \\( m \\) as \\( \mathbb{Z} / m\mathbb{Z} \\). We say that \\( a + m\mathbb{Z} \\) is invertible in \\( \mathbb{Z} / m\mathbb{Z} \\) if and only if there exists a solution for \\( ax \equiv 1 \pmod{m} \\).*

---

**Example:** In \\( \mathbb{Z}/5\mathbb{Z} \\), \\( 3 + 5\mathbb{Z} \\) is invertible because \\( \gcd(3, 5) = 1 \\) (since \\( 3\cdot2 \equiv 1 \pmod{5} \\)). However, in \\( \mathbb{Z}/6\mathbb{Z} \\), \\( 3 + 6\mathbb{Z} \\) is not invertible because \\( \gcd(3, 6) = 3 \neq 1 \\).

### Lemma 2.2.1 

---

*If \\( a \\) and \\( b \\) are coprime, the residues of \\( a \\), \\( 2a \\), \\( 3a \\), ..., \\( (b-1)a \\) modulo \\( b \\) are all distinct.*

---

**Proof:** Suppose, for contradiction, that there exist \\( x, y \in \\{1, 2, \dots, b-1\\} \\) with \\( x \neq y \\) such that \\( xa \equiv ya \pmod{b} \\). Then \\( (x-y)a \equiv 0 \pmod{b} \\), implying \\( b \mid (x-y)a \\). Since \\( a \\) and \\( b \\) are coprime, we must have \\( b \mid (x-y) \\). However, \\( |x-y| < b \\), so this is only possible if \\( x = y \\), contradicting our assumption. Therefore, all residues must be distinct.

### Theorem 2.2.2

---

*For any integers \\( a \\) and \\( b \\), the equation \\( ax + by = 1 \\) has a solution in integers \\( x \\) and \\( y \\) if and only if \\( a \\) and \\( b \\) are coprime.*

---

**Proof:**  

- **(\\( \Rightarrow \\))** Suppose \\( a \\) and \\( b \\) are not coprime, let \\( d = \gcd(a,b) > 1 \\). Then \\( d \mid a \\) and \\( d \mid b \\), so \\( d \mid (ax + by) \\) for any integers \\( x \\) and \\( y \\). Thus, \\( ax + by \neq 1 \\) for any \\( x \\) and \\( y \\).
- **(\\( \Leftarrow \\))** Suppose \\( a \\) and \\( b \\) are coprime. By the previous lemma, the residues of \\( a \\), \\( 2a \\), ..., \\( (b-1)a \\) modulo \\( b \\) are all distinct. Therefore, there exists an \\( m \in \\{1, 2, ..., b-1\\} \\) such that \\( ma \equiv 1 \pmod{b} \\). This means there exists an integer \\( n \\) such that \\( ma = bn + 1 \\), or equivalently, \\( ma - bn = 1 \\).

### Theorem: Bézout's Identity

---

*For any integers \\( a \\) and \\( b \\), we have:
\\(
a\mathbb{Z} + b \mathbb{Z} = \gcd(a, b)\mathbb{Z}
\\)*

---

**Proof:** 

This statement is equivalent to proving that \\( ax + by = c \\) has an integer solution if and only if \\( c \\) is a multiple of \\( \gcd(a,b) \\).

- **(\\( \Rightarrow \\))** If \\( ax + by = c \\) for some integers \\( x \\) and \\( y \\), then \\( \gcd(a,b) \mid a \\) and \\( \gcd(a,b) \mid b \\), so \\( \gcd(a,b) \mid (ax + by) = c \\).
- **(\\( \Leftarrow \\))** Let \\( c = k\gcd(a,b) \\) for some integer \\( k \\). We can write \\( a = p\gcd(a,b) \\) and \\( b = q\gcd(a,b) \\), where \\( p \\) and \\( q \\) are coprime. By the previous theorem, there exist integers \\( m \\) and \\( n \\) such that \\( pm + qn = 1 \\). Multiplying both sides by \\( k\gcd(a,b) \\), we get:
  \\[
  akm + bkn = c
  \\]
  Thus, \\( x = km \\) and \\( y = kn \\) are integer solutions to \\( ax + by = c \\).

This theorem implies that for any integers \\( a \\), \\( b \\), and \\( n \\), the equation \\( ax + by = n \\) has an integer solution if and only if \\( \gcd(a,b) \mid n \\).

### Theorem 2.2.3

---

*\\( a + m\mathbb{Z} \\) is invertible in \\( \mathbb{Z}/m\mathbb{Z} \\) if and only if \\( \gcd(a,m) = 1 \\).*

---

**Proof:**  

- **(\\( \Rightarrow \\))** Suppose \\( a + m\mathbb{Z} \\) is invertible in \\( \mathbb{Z}/m\mathbb{Z} \\). Then there exists an integer \\( x \\) such that \\( ax \equiv 1 \pmod{m} \\). Let \\( g = \gcd(a,m) \\). Since \\( g \mid ax \\) and \\( g \mid (ax - 1) \\), we must have \\( g \mid 1 \\), so \\( g = 1 \\).
- **(\\( \Leftarrow \\))** Suppose \\( \gcd(a,m) = 1 \\). By Bézout's identity, there exist integers \\( x \\) and \\( y \\) such that \\( ax + my = 1 \\). Thus, \\( x + m\mathbb{Z} \\) is the multiplicative inverse of \\( a + m\mathbb{Z} \\) in \\( \mathbb{Z}/m\mathbb{Z} \\).

### Definition: Residue Class Ring

---

*\\( (\mathbb{Z} / m \mathbb{Z}, +, \cdot) \\) is a commutative ring where \\( 1 + m \mathbb{Z} \\) is the multiplicative identity element. This ring is called the residue class ring modulo \\( m \\).*

---

**Example:** \\( \mathbb{Z}/4\mathbb{Z} = \\{0 + 4\mathbb{Z}, 1 + 4\mathbb{Z}, 2 + 4\mathbb{Z}, 3 + 4\mathbb{Z}\\} \\).

### Definition: Primitive Residue Class

---

*A residue class \\( a + m\mathbb{Z} \\) is called primitive if \\( \gcd(a, m) = 1 \\).*

---

**Example:** In \\( \mathbb{Z}/6\mathbb{Z} \\), the primitive residue classes are \\( 1 + 6\mathbb{Z} \\) and \\( 5 + 6\mathbb{Z} \\).

### Theorem 2.2.4

---

*A residue ring \\( \mathbb{Z} / m\mathbb{Z} \\) is a field if and only if \\( m \\) is a prime number.*

---

**Proof:** TBD

**Example:** For \\( m = 5 \\), \\( \mathbb{Z}/5\mathbb{Z} = \\{0 + 5\mathbb{Z}, 1 + 5\mathbb{Z}, 2 + 5\mathbb{Z}, 3 + 5\mathbb{Z}, 4 + 5\mathbb{Z}\\} \\) forms a field because 5 is a prime number, and every non-zero element has a multiplicative inverse.

### Definition: Primitive Residue Class Group

---

*The group of all primitive residue classes modulo \\( m \\) is called the primitive residue class group, denoted by \\( (\mathbb{Z}/m\mathbb{Z})^{\times} \\).*

---

**Example:** For \\( m = 8 \\), the set of all primitive residue classes is \\( (\mathbb{Z}/8\mathbb{Z})^{\times} = \\{1 + 8\mathbb{Z}, 3 + 8\mathbb{Z}, 5 + 8\mathbb{Z}, 7 + 8\mathbb{Z}\\} \\). These are the integers less than 8 that are coprime to 8 (i.e., \\(\gcd(a, 8) = 1\\)).

Contrast this with \\( m = 9 \\). The primitive residue class group is \\((\mathbb{Z}/9\mathbb{Z})^{\times} = \\{1 + 9\mathbb{Z}, 2 + 9\mathbb{Z}, 4 + 9\mathbb{Z}, 5 + 9\mathbb{Z}, 7 + 9\mathbb{Z}, 8 + 9\mathbb{Z}\\}\\), as these are the integers less than 9 that are coprime to 9.

### Definition: Euler's Totient Function

---

*Euler's totient function \\(\phi(m)\\) is equal to the order of the primitive residue class group modulo \\(m\\), which is the number of integers less than \\(m\\) and coprime to \\(m\\).*

---

**Example:** For \\(m = 12\\), \\(\phi(12) = 4\\) because there are 4 integers less than 12 that are coprime to 12: \\(\{1, 5, 7, 11\}\\).

For \\(m = 10\\), \\(\phi(10) = 4\\), as there are also 4 integers less than 10 that are coprime to 10: \\(\{1, 3, 7, 9\}\\).

### Definition: Order of an Element within a Group

---

*The order of \\(g \in G\\) is the smallest natural number \\(e\\) such that \\(g^{e} = 1\\). We denote it as \\(\text{order}_g(g)\\)or \\(\text{order } g\\).*

---

**Example:** In \\((\mathbb{Z}/7\mathbb{Z})^{\times}\\), the element 3 has order 6 because \\(3^6 \bmod 7 = 1\\). In other words, \\(3 \times 3 \times 3 \times 3 \times 3 \times 3 \bmod 7 = 1\\), and 6 is the smallest such exponent.

### Definition: Subgroup

---

*A subset \\(U \subseteq G\\) is a subgroup of \\(G\\) if \\(U\\) itself is a group under the operation of \\(G\\).*

---

**Example:** Consider \\((\mathbb{Z}/8\mathbb{Z})^{\times} = \\{1 + 8\mathbb{Z}, 3+ 8\mathbb{Z}, 5+ 8\mathbb{Z}, 7+ 8\mathbb{Z}\\}\\) under multiplication modulo 8. The subset \\(\{1+ 8\mathbb{Z}, 7+ 8\mathbb{Z}\}\\) forms a subgroup because it satisfies the group properties: closed under multiplication, contains the identity element (1), and every element has an inverse (\\(7 \times 7 \equiv 1 \bmod 8\\)).

### Definition: Subgroup Generated by \\(g\\)

---

*The set \\(\\{g^{k} : k \in \mathbb{Z}\\}\\), for some element \\(g \in G\\), forms a subgroup of \\(G\\) and is called the subgroup generated by \\(g\\), denoted by \\(\langle g \rangle\\).*

---

**Example:** Consider the group \\((\mathbb{Z}/7\mathbb{Z})^{\times} = \\{1+ 7\mathbb{Z}, 2+ 7\mathbb{Z}, 3+ 7\mathbb{Z}, 4+ 7\mathbb{Z}, 5+ 7\mathbb{Z}, 6+ 7\mathbb{Z}\\}\\) under multiplication modulo 7. If we take \\(g = 3\\), then \\(\langle 3 +7\mathbb{Z} \rangle =\\) \\(\\{3^1+7\mathbb{Z}, 3^2+7\mathbb{Z}, 3^3+7\mathbb{Z}, 3^4+7\mathbb{Z}, 3^5+7\mathbb{Z}, 3^6+7\mathbb{Z}\\} \bmod 7 =\\) \\(\\{3+7\mathbb{Z}, 2+7\mathbb{Z}, 6+7\mathbb{Z}, 4+7\mathbb{Z}, 5+7\mathbb{Z}, 1+7\mathbb{Z}\\}\\), which forms a subgroup generated by 3. This subgroup contains all elements of \\((\mathbb{Z}/7\mathbb{Z})^{\times}\\), making 3 a generator of the entire group.

If \\(g\\) has a finite order \\(e\\), we have that \\(\langle g \rangle = \\{g^{k}: 0 \leq k \leq e\\}\\), meaning \\(e\\) is the order of \\(\langle g \rangle\\).

### Definition: Cyclic Group

---

*A group \\(G\\) is called a cyclic group if there exists an element \\(g \in G\\) such that \\(G = \langle g \rangle\\). In this case, \\(g\\) is called a generator of \\(G\\).*

---

**Example:** The group \\((\mathbb{Z}/6\mathbb{Z})^{\times} = \\{1+6\mathbb{Z}, 5+6\mathbb{Z}\\}\\) under multiplication modulo 6 is a cyclic group. In this case, both 1 and 5 are generators of the group because \\(\langle 5 +6\mathbb{Z} \rangle = \\{(5^1 \bmod 6)+6\mathbb{Z} = 5 +6\mathbb{Z}, (5^2 \bmod 6)+6\mathbb{Z} = 1+6\mathbb{Z}\\}\\). Since 5 generates all the elements of the group, \\(G\\) is cyclic.

### Theorem 2.2.5

---

*If \\(G\\) is a finite cyclic group, it has \\(\phi(|G|)\\)generators, and each generator has order \\(|G|\\).*

---

**Proof**: TBD

**Example:** Consider the group \\((\mathbb{Z}/8\mathbb{Z})^{\times} = \\{1+8\mathbb{Z}, 3+8\mathbb{Z}, 5+8\mathbb{Z}, 7+8\mathbb{Z}\\}\\). This group is cyclic, and \\(\phi(8) = 4\\). The generators of this group are \\(\\{1+8\mathbb{Z}, 3+8\mathbb{Z}, 5+8\mathbb{Z}, 7+8\mathbb{Z}\\}\\), each of which generates the entire group when raised to successive powers modulo 8. Each generator has the same order, which is \\(|G| = 4\\).

### Theorem 2.2.6

---

*If \\(G\\) is a finite cyclic group, the order of any subgroup of \\(G\\) divides the order of \\(G\\).*

---

**Proof**: TBD

**Example:** Consider the cyclic group \\((\mathbb{Z}/6\mathbb{Z})^{\times} = \\{1+6\mathbb{Z}, 5+6\mathbb{Z}\\}\\) under multiplication modulo 6. If we take the subgroup \\(\langle 5+6\mathbb{Z} \rangle = \\{1+6\mathbb{Z}, 5+6\mathbb{Z}\\}\\), this is a subgroup of order 2, and 2 divides the order of the original group, which is 6. This theorem generalizes this property: for any subgroup of a cyclic group, its order divides the order of the group.

### Theorem: Fermat's Little Theorem

---

*If \\(\gcd(a, m) = 1\\), then \\(a^{\phi(m)} \equiv 1 \pmod{m}\\).*

---

**Proof**: TBD

**Example:** Take \\(a = 2\\) and \\(m = 5\\). Since \\(\gcd(2, 5) = 1\\), Fermat's Little Theorem tells us that \\(2^{\phi(5)} = 2^4 \equiv 1 \bmod 5\\). Indeed, \\(2^4 = 16\\) and \\(16 \bmod 5 = 1\\).

This theorem suggests that \\(a^{\phi(m) - 1} + m \mathbb{Z}\\) is the inverse residue class of \\(a + m \mathbb{Z}\\).

### Theorem 2.2.7

---

*The order of any element in a group divides the order of the group.*

---

**Proof**: TBD

**Example:** In the group \\((\mathbb{Z}/7\mathbb{Z})^{\times}\\), consider the element \\(3 + 7\mathbb{Z}\\). The order of \\(3 + 7\mathbb{Z}\\) is 6, as \\(3^6 \equiv 1 \bmod 7\\). The order of the group itself is also 6, and indeed, the order of the element divides the order of the group.

### Theorem: Generalization of Fermat's Little Theorem

---

*For any element \\(g \in G\\), we have \\(g^{|G|} = 1\\).*

---

**Proof**: TBD

**Example:** In the group \\((\mathbb{Z}/7\mathbb{Z})^{\times}\\), for any element \\(g\\), such as \\(g = 3 + 7\mathbb{Z}\\), we have \\(3^6 \equiv 1 \bmod 7\\). This holds for any \\(g \in (\mathbb{Z}/7\mathbb{Z})^{\times}\\) because the order of the group is 6. Thus, \\(g^{|G|} = 1\\) is satisfied.




