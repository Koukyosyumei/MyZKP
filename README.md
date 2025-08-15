```
███╗   ███╗  ██╗   ██╗  ███████╗  ██╗  ██╗  ██████╗   🦀
████╗ ████║  ╚██╗ ██╔╝  ╚══███╔╝  ██║ ██╔╝  ██╔══██╗ 🦀
██╔████╔██║   ╚████╔╝     ███╔╝   █████╔╝   ██████╔╝ 
██║╚██╔╝██║    ╚██╔╝     ███╔╝    ██╔═██╗   ██╔═══╝    🦀
██║ ╚═╝ ██║     ██║     ███████╗  ██║  ██╗  ██║    🦀    
╚═╝     ╚═╝     ╚═╝     ╚══════╝  ╚═╝  ╚═╝  ╚═╝      🦀
```

# 🚀 Building Zero Knowledge Proofs from Scratch in Rust

**MyZKP** is a Rust implementation of zero-knowledge protocols built entirely from scratch! This project serves as an educational resource for understanding and working with zero-knowledge proofs. 

- 🌟 Explore the accompanying textbook: [**Book of MyZKP**](https://koukyosyumei.github.io/MyZKP/).
- 📢 Contributions and feedback are encouraged to help make this project even better!


> [!WARNING]  
> This repository is a work in progress and may contain bugs or inaccuracies. Contributions and feedback are welcome!

## 📚 About MyZKP

MyZKP is a growing library that provides:

- A step-by-step guide to the theoretical foundations of ZKPs, including number theory, elliptic curves, and field arithmetic.
- Implementation of core primitives for ZKP protocols.
- A solid base for developers and researchers to learn, experiment, or build their own ZKP-based systems.
  
💡 Whether you're a cryptography enthusiast, a Rustacean, or a student, MyZKP is for you!


## 📖 Educational Modules

**🧮 Basic of Number Theory**

- 📝 [Computation Rule and Properties](https://koukyosyumei.github.io/MyZKP/number_theory/subsec1.html)
- ⚙️ [Semigroup, Group, Ring, and Field](https://koukyosyumei.github.io/MyZKP/number_theory/subsec2.html)
- 🔢 [Polynomials](https://koukyosyumei.github.io/MyZKP/number_theory/subsec3.html)
- 🌐 [Galois Field](https://koukyosyumei.github.io/MyZKP/number_theory/subsec4.html)
- 📈 [Elliptic Curve](https://koukyosyumei.github.io/MyZKP/number_theory/subsec5.html)
- 🔗 [Pairing](https://koukyosyumei.github.io/MyZKP/number_theory/subsec6.html)
- 🤔 [Useful Assumptions](https://koukyosyumei.github.io/MyZKP/number_theory/subsec7.html)

**🔒 Basic of zk-SNARKs**

- ⚡ [Arithmetization](https://koukyosyumei.github.io/MyZKP/zksnark/subsec1.html)
- 🛠️ [Proving Single Polynomial](https://koukyosyumei.github.io/MyZKP/zksnark/subsec2.html)
- 🐍 [Bringing It All Together: SNARK](https://koukyosyumei.github.io/MyZKP/zksnark/subsec3.html)
  
**🌟 Basic of zk-STARKs**

- 🌈 [FRI](https://koukyosyumei.github.io/MyZKP/zkstark/subsec1.html)

**💻 Basic of zkVM**

- ✍️ TBD

## 🛠️ Code Reference

|Module       | Submodule   | Description  |📂 Path                                      |
|-------------|-------------|--------------|-------------------------------------------------- |
| **albebra** | `ring`      | Defines Ring  | [ring.rs](./myzkp/src/modules/algebra/ring.rs)           |
|             | `field`     | Defines Finite Field | [field.rs](./myzkp/src/modules/algebra/field.rs)           |
|             | `efield`    | Field Extension (Galois Field) |[efield.rs](./myzkp/src/modules/algebra/efield.rs)           |
|             | `Polynomial`| Polynomial Arithmetic | [polynomial.rs](./myzkp/src/modules/algebra/polynomial.rs)|
|             | `curve`     | Elliptic curve operations | [curve](./myzkp/src/modules/algebra/curve/)           |
| **Arithmetization** | `r1cs` | Rank One Constraint System | [r1cs.rs](./myzkp/src/modules/arithmetization/r1cs.rs) |
|             | `qap` | Quadratic Arithmetic Program | [qap.rs](./myzkp/src/modules/arithmetization/qap.rs) |
| **Polynomial Commitment** | `kzg`| [Kate, Aniket, Gregory M. Zaverucha, and Ian Goldberg. "Constant-size commitments to polynomials and their applications." International conference on the theory and application of cryptology and information security. Berlin, Heidelberg: Springer Berlin Heidelberg, 2010.](https://link.springer.com/chapter/10.1007/978-3-642-17373-8_11) | [kzg.rs](./myzkp/src/modules/algebra/kzg.rs) |
|                           | `gemini`| [Bootle, Jonathan, et al. "Gemini: Elastic SNARKs for diverse environments." Annual International Conference on the Theory and Applications of Cryptographic Techniques. Cham: Springer International Publishing, 2022.](https://link.springer.com/chapter/10.1007/978-3-031-07085-3_15) | [gemini.rs](./myzkp/src/modules/algebra/gemini.rs) |
| **zkSNARKs**| `tutorial_single_polynomial` | [Petkus, Maksym. "Why and how zk-snark works." arXiv preprint arXiv:1906.07221 (2019).](https://arxiv.org/abs/1906.07221) | [tutorial_single_polynomial](./myzkp/src/modules/zksnark/tutorial_single_polynomial/)                                   |
|             | `tutorial_snark` | [Petkus, Maksym. "Why and how zk-snark works." arXiv preprint arXiv:1906.07221 (2019).](https://arxiv.org/abs/1906.07221) | [tutorial_snark](./myzkp/src/modules/zksnark/tutorial_snark/) |
|             | `pinocchio` | [Parno, Bryan, et al. "Pinocchio: Nearly practical verifiable computation." Communications of the ACM 59.2 (2016): 103-112.](https://dl.acm.org/doi/abs/10.1145/2856449) | [pinocchio.rs](./myzkp/src/modules/zksnark/pinocchio.rs) |
| **zkSTARKs**| `fri`       | [Ben-Sasson, Eli, et al. "DEEP-FRI: sampling outside the box improves soundness." arXiv preprint arXiv:1903.12243 (2019).](https://arxiv.org/abs/1903.12243) | [fri.rs](./myzkp/src/modules/zkstark/fri.rs) |

## 🤝 Contributions are Welcome!

We welcome your ideas, feedback, and contributions!

💡 Here are ways to contribute:

1. *Report Issues*: Found a bug or have a suggestion? Open an issue!
2. *Submit Pull Requests*: Have code improvements? Feel free to submit them.
3. *Write Documentation*: Help us make the educational modules clearer.