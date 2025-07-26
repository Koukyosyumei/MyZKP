# Basics of zk-STARK

zk-Stark stands for **S**calable **T**ransparent **AR**gument of **K**nowledge. Its key properties are:

- Scalability: Prover time grows as \\(\mathcal{O}(T\log{T})\\) and verifier time as \\(\mathcal{O}(\log{T})\\), where \\(T\\) is the cost of executing the statement being proven.
- Transparency: Unlike many SNARK systems, zk‑STARKs require no trusted setup.
- Post‑Quantum Security: Security depends solely on the collision and preimage resistance of cryptographic hash functions, making zk‑STARKs resilient against quantum‑computing attacks.

In this tutorial, we will unpack and implement the fundamental building blocks of zk‑STARKs. Our presentation draws heavily on the [Anatomy of a STARK](https://aszepieniec.github.io/stark-anatomy/) guide by Alan Szepieniec and the [STARK 101](https://starkware.co/stark-101/) series from StarkWare.

