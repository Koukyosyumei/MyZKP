use num_traits::{One, Zero};
use std::ops::{Div, Mul, Rem, Sub};

pub fn extended_euclidean<F>(a: F, b: F) -> (F, F, F)
where
    F: Clone + PartialEq + Sub<Output = F> + Div<Output = F> + Rem<Output = F> + Zero + One,
    for<'a> &'a F: Sub<&'a F, Output = F>
        + Div<&'a F, Output = F>
        + Rem<&'a F, Output = F>
        + Mul<&'a F, Output = F>,
{
    let mut r0 = a;
    let mut r1 = b;
    let mut s0 = F::one();
    let mut s1 = F::zero();
    let mut t0 = F::zero();
    let mut t1 = F::one();

    while !r1.is_zero() {
        let q = &r0 / &r1;
        let r = &r0 % &r1;
        r0 = r1;
        r1 = r;
        let new_s = s0 - &q * &s1;
        s0 = s1;
        s1 = new_s;
        let new_t = t0.clone() - &q * &t1;
        t0 = t1;
        t1 = new_t;
    }

    (r0, s0, t0)
}
