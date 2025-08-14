# Polynomial Commitment

Good point — let’s start with a short, formal definition of *polynomial commitments* and then embed that as the intro to the KZG write-up.

# Polynomial commitment (definition)

A **polynomial commitment scheme** lets a prover commit to a polynomial \\(f(X)\\) so that later the prover can (1) reveal the value \\(f(u)\\) at any point \\(u\\) together with a *short* proof (witness), and (2) the verifier can efficiently check that the revealed value is correct relative to the original commitment. A scheme typically provides four algorithms:

* **Setup**\\((1^\lambda, D)\\): (possibly trusted) produces public parameters \\(\text{PK}\\) and possibly a secret key \\(\text{SK}\\). \\(D\\) is a public upper bound on polynomial degree.
* **Commit**\\((\text{PK}, f)\to C\\): the prover computes a short commitment \\(C\\) to polynomial \\(f\\).
* **Open**\\((\text{PK}, f, u)\to (y, W)\\): the prover computes \\(y=f(u)\\) and a witness \\(W\\) that proves \\(y\\) is the value of \\(f\\) at \\(u\\).
* **Verify**\\((\text{PK}, C, u, y, W)\to\\{\text{accept},\text{reject}\\}\\): the verifier checks the witness and either accepts or rejects.

**Correctness.** For honestly generated keys and honest prover, if \\(C\leftarrow\text{Commit}(f)\\) and \\((y,W)\leftarrow\text{Open}(f,u)\\), then \\(\text{Verify}(\text{PK},C,u,y,W)\\) must accept.

**Security properties (informal).**

* **Binding:** After publishing \\(C\\), the committer cannot (except with negligible probability under standard assumptions) produce two different openings \\((u,y,W)\\) and \\((u,y',W')\\) with \\(y\ne y'\\) that both verify. This prevents equivocation.
* **Hiding (optional):** Some schemes hide the polynomial (or its evaluations) from the verifier; others are not hiding by default. Hiding can be statistical or computational depending on the scheme and added randomness.
* **Succinctness / efficiency:** Commitments and opening proofs should be short (ideally constant size independent of polynomial degree), and verification should be much cheaper than re-evaluating \\(f\\) from its coefficients.
