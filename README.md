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

- ⚡ [Arithmetization](https://koukyosyumei.github.io/MyZKP/zksnark/subsec2.html)
- 🛠️ [Proving Single Polynomial](https://koukyosyumei.github.io/MyZKP/zksnark/subsec3.html)
- 🐍 [Bringing It All Together: SNARK](https://koukyosyumei.github.io/MyZKP/zksnark/subsec4.html)
  
**🌟 Basic of zk-STARKs**

- ✍️ TBD

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
| **zkSNARKs**| `tutoriak_single_polynomial` | | [tutorial_single_polynomial](./myzkp/src/modules/zksnark/tutorial_single_polynomial/)                                   |
|             | `tutorial_snark` | | [tutorial_snark](./myzkp/src/modules/zksnark/tutorial_snark/) |

## 🤝 Contributions are Welcome!

We welcome your ideas, feedback, and contributions!

💡 Here are ways to contribute:

1. *Report Issues*: Found a bug or have a suggestion? Open an issue!
2. *Submit Pull Requests*: Have code improvements? Feel free to submit them.
3. *Write Documentation*: Help us make the educational modules clearer.