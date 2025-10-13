use num_traits::identities::One;
use myzkp::modules::algebra::field::{FiniteFieldElement, ModEIP197};

type F = FiniteFieldElement<ModEIP197>;

fn main() {
    let a = F::one();
}
