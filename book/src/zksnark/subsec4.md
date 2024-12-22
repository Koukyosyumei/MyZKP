# Bringing It All Together

We first recap the previous sections. First the relationship between inputs and outputs of any prgram can be represented as a rank-one constraint system (R1CS) as follows:

\\[
  (L \cdot a) \circ (R \cdot a) - O \cdot a = 0  
\\]

, where \\(a\\) is the concatenation of all inputs, outputs, and intermediate values (witness). Then, the statement, "I know the input values \\(x\\) that make the program returns the output values \\(y\\)", can be converted into "I know \\(a\\), whose outputs part is \\(y\\), that satisfies the constraint system corresponding to the program". 

Although separately checking each constraint, corresponding to each row in the above R1CS, is not efficient, we can convert those vector-equivalence test into the polynomial-equivalence test.

## First Protocol Idea

The simplest protocol based on the previous chapter is as follos:

**Protocol (Setup)**

- **Interpolated Polynomial:** Construct \\(\\{\ell_i, r_i, o_i\\}_{i\in[d]}\\) from \\(L\\), \\(R\\), and \\(O\\), respectively.
- **Target Polynomial:** \\(t(x) = (x-1)(x-2) \cdots (x-m)\\)
- **Secret Seed:** A trusted party generates the random value \\(s\\) and \\(\alpha\\).
- **Proof Key:** Provided to the prover
  - \\(\\{g^{\ell_i(s)},g^{r_i(s)},g^{o_i(s)}\\}_{i\in[d]}\\)
  - \\(\\{g^{\alpha \ell_i(s)},g^{\alpha r_i(s)},g^{\alpha o_i(s)}\\}_{i\in[d]}\\)
  - \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Verification Key:**
  - \\(g^{t(s)}\\)
  - \\(g^{\alpha}\\)
- After distribution, the original \\(s\\) and \\(\alpha\\) values are securely destroyed.

**Protocol (Proving)**

- Execute the program and get the witness vector \\(w\\).
- Compute the linear-combinations of polynomials
  - \\(\ell(x) = \sum_{i=1}^{d} w_i \ell_{i}(x)\\)
  - \\(r(x) = \sum_{i=1}^{d} w_i r_{i}(x)\\)
  - \\(o(x) = \sum_{i=1}^{d} w_i o_{i}(x)\\)
- Compute \\(h(x) = \frac{\ell(x) r(x) - o(x)}{t(x)}\\)
- Evaluate each polynomial at \\(s\\).
  - \\(g^{\ell(s)} = \prod^{d}_{i=1} (g^{\ell_i(s)})^{w_i} \\)
  - \\(g^{r(s)} = \prod^{d}_{i=1} (g^{r_i(s)})^{w_i} \\)
  - \\(g^{o(s)} = \prod^{d}_{i=1} (g^{o_i(s)})^{w_i} \\)
- Evaluate each shifted polynomial at \\(s\\).
  - \\(g^{\alpha \ell(s)} = \prod^{d}_{i=1} (g^{\alpha \ell_i(s)})^{w_i} \\)
  - \\(g^{\alpha r(s)} = \prod^{d}_{i=1} (g^{\alpha r_i(s)})^{w_i} \\)
  - \\(g^{\alpha o(s)} = \prod^{d}_{i=1} (g^{\alpha o_i(s)})^{w_i} \\)
- Calculate \\(g^{h(s)}\\) from \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- Proof: \\((g^{\ell(s)}, g^{r(s)}, g^{o(s)}, g^{\alpha \ell(s)}, g^{\alpha r(s)}, g^{\alpha o(s)}, g^{h(s)})\\)

**Protocol (Verification)**

- Parse proof as \\((g^{\ell}, g^r, g^o, g^{\ell'}, g^{r'}, g^{o'}, g^{h})\\)
- Check polynomial restrictions
  - \\(e(g^{\ell}, g^{\alpha}) = e(g^{\ell'}, g)\\)
  - \\(e(g^{r}, g^{\alpha}) = e(g^{r'}, g)\\)
  - \\(e(g^{o}, g^{\alpha}) = e(g^{o'}, g)\\)
- Validity check
  - \\(e(g^{\ell}, g^{r}) = e(g^t,g^h) \cdot e(g^o, g)\\)

**Vulnerability**

One ciriticl problem of the above protocol is that its check may pass even if \\(g^{\ell}\\) is computed from not \\(\\{g^{\ell_i(s)}\\}_{i \in [d]}\\) but either \\(\\{g^{r_i(s)}\\} _{i \in [d]}\\) or \\(\\{g^{o_i(s)}\\} _{i \in [d]}\\), or their combination. The same phenomenon is also true for \\(g^{r}\\) and \\(g^{o}\\). For example, if the prover sends \\((g^{r(s)}, g^{\ell(s)}, g^{o(s)}, g^{\alpha r(s)}, g^{\alpha \ell(s)}, g^{\alpha o(s)}, g^{h(s)})\\) insteads of \\((g^{\ell(s)}, g^{r(s)}, g^{o(s)}, g^{\alpha \ell(s)}, g^{\alpha r(s)}, g^{\alpha o(s)}, g^{h(s)})\\), all of the verification checks still pass, although the proved statement is different from the original one.

## Second Protocol: Non-Interchangibility

To solve the above interchangeability problem, we should use the different \\(\alpha\\)-shift for \\(\ell\\), \\(r\\), and \\(o\\).

**Protocol (Setup)**

- **Interpolated Polynomial:** Construct \\(\\{\ell_i, r_i, o_i\\}_{i\in[d]}\\) from \\(L\\), \\(R\\), and \\(O\\), respectively.
- **Target Polynomial:** \\(t(x) = (x-1)(x-2) \cdots (x-m)\\)
- **Secret Seed:** A trusted party generates the random value \\(s\\), *\\(\alpha_{\ell}\\), \\(\alpha_r\\), and \\(\alpha_o\\)*.
- **Proof Key:** Provided to the prover
  - \\(\\{g^{\ell_i(s)},g^{r_i(s)},g^{o_i(s)}\\}_{i\in[d]}\\)
  - *\\(\\{g^{\alpha_{\ell} \ell_i(s)},g^{\alpha_{r} r_i(s)},g^{\alpha_{o} o_i(s)}\\}_{i\in[d]}\\)*
  - \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Verification Key:**
  - \\(g^{t(s)}\\)
  - *\\(g^{\alpha_{\ell}}\\), \\(g^{\alpha_{r}}\\), \\(g^{\alpha_{o}}\\)*
- After distribution, the original \\(s\\), \\(\alpha_{\ell}\\), \\(\alpha_r\\), and \\(\alpha_o\\) values are securely destroyed.

**Protocol (Proving)**

- Execute the program and get the witness vector \\(w\\).
- Compute the linear-combinations of polynomials
  - \\(\ell(x) = \sum_{i=1}^{d} w_i \ell_{i}(x)\\)
  - \\(r(x) = \sum_{i=1}^{d} w_i r_{i}(x)\\)
  - \\(o(x) = \sum_{i=1}^{d} w_i o_{i}(x)\\)
- Compute \\(h(x) = \frac{\ell(x) r(x) - o(x)}{t(x)}\\)
- Evaluate each polynomial at \\(s\\).
  - \\(g^{\ell(s)} = \prod^{d}_{i=1} (g^{\ell_i(s)})^{w_i} \\)
  - \\(g^{r(s)} = \prod^{d}_{i=1} (g^{r_i(s)})^{w_i} \\)
  - \\(g^{o(s)} = \prod^{d}_{i=1} (g^{o_i(s)})^{w_i} \\)
- Evaluate each shifted polynomial at \\(s\\).
  - *\\(g^{\alpha_{\ell} \ell(s)} = \prod^{d}_{i=1} (g^{\alpha _{\ell} \ell_i(s)})^{w_i} \\)*
  - *\\(g^{\alpha_{r} r(s)} = \prod^{d}_{i=1} (g^{\alpha _{r} r_i(s)})^{w_i} \\)*
  - *\\(g^{\alpha_{o} o(s)} = \prod^{d}_{i=1} (g^{\alpha _{o} o_i(s)})^{w_i} \\)*
- Calculate \\(g^{h(s)}\\) from \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- Proof: \\((g^{\ell(s)}, g^{r(s)}, g^{o(s)},\\) *\\(g^{\alpha_{\ell} \ell(s)}, g^{\alpha_{r} r(s)}, g^{\alpha_{o} o(s)},\\)* \\(g^{h(s)})\\)

**Protocol (Verification)**

- Parse proof as \\((g^{\ell}, g^r, g^o, g^{\ell'}, g^{r'}, g^{o'}, g^{h})\\)
- Check polynomial restrictions
  - *\\(e(g^{\ell}, g^{\alpha_{\ell}}) = e(g^{\ell'}, g)\\)*
  - *\\(e(g^{r}, g^{\alpha_{r}}) = e(g^{r'}, g)\\)*
  - *\\(e(g^{o}, g^{\alpha_{o}}) = e(g^{o'}, g)\\)*
- Validity check
  - \\(e(g^{\ell}, g^{r}) = e(g^t,g^h) \cdot e(g^o, g)\\)

**Vulnerability**

The current protocol still does not force the consistency of each variable. In other words, each variable \\(w_i\\) can be different values in each \\(\ell\\), \\(r\\), and \\(o\\), since the verification check is done separetely.

## Thrid Protocol: Variable Consistency

To achive the variable-consistency, we use a kind of checksum. Specifically, we first draw a new random value \\(\beta\\) and defines a checksum of \\(w_i\\) as \\(g^{\beta(\ell_{i}(s) + r_i(s) + o_i(s))}\\). Let \\(w _{\ell, i}\\), \\(w _{r,i}\\), \\(w _{o,i}\\), and \\(w _{\beta,i}\\) be the \\(i\\)-th value of the witness vector for \\(\ell\\), \\(r\\), \\(o\\) and the checksum, respectively. Then, if all of them are true, the following equation holds:

\\[
e(g^{w _{\ell, i} \ell_i(s)} g^{w _{r, i} r_i(s)} g^{w _{o, i} o_i(s)}, g^{\beta}) = e(g^{w _{\beta, i} \beta(\ell _{i}(s) + r _{i}(s) + o _{i}(s))}, g)
\\]

Unfortunately, this condition is not equivalent. For example, suppose \\(\ell_i(x) = r_i(x)\\). Then, we have the following:

\begin{align*}
  &\beta(w _{\ell, i} \ell_i(s) + w _{r, i} r_i(s) + w _{o, i} o_i(s)) = w _{\beta, i} \beta (\ell _{i}(s) + r _{i}(s) + o _{i}(s)) \\\\
  \iff &\beta(w _{\ell, i} \ell_i(s) + w _{r, i} \ell_i(s) + w _{o, i} o_i(s)) = w _{\beta, i} \beta (2\ell _{i}(s) + o _{i}(s))
\end{align*}

Then, for arbitrary \\(w _{r,i}\\) and \\(w _{o,i}\\), the above equation holds if we set \\(w _{\beta, i} = w _{o, i}\\) and \\(w _{\ell, i} = 2 w _{o, i} - w _{r, i}\\).

To surrogate this problem, we have to use different \\(\beta\\) for each of \\(\ell\\), \\(r\\) and \\(o\\). THen, we need to check the following equation:

\\[
e(g^{w _{\ell, i} \ell _{i}(s)}, g^{\beta _{\ell}}) \cdot e(g^{w _{r, i} r _{i}(s)}, g^{\beta _{r}}) \cdot e(g^{w _{o, i} o _{i}(s)}, g^{\beta _{o}}) = e(g^{w _{\beta, i}(\beta _{\ell} \ell _{i}(s) + \beta _{r} r _{i}(s) + \beta _{o} o _{i}(s))}, g)   
\\]

The new protocol using the above variable-consistency check is as follows:

**Protocol (Setup)**

- **Interpolated Polynomial:** Construct \\(\\{\ell_i, r_i, o_i\\}_{i\in[d]}\\) from \\(L\\), \\(R\\), and \\(O\\), respectively.
- **Target Polynomial:** \\(t(x) = (x-1)(x-2) \cdots (x-m)\\)
- **Secret Seed:** A trusted party generates the random value \\(s\\), \\(\alpha_{\ell}\\), \\(\alpha_r\\), \\(\alpha_o\\), *\\(\beta_{\ell}\\), \\(\beta_{r}\\), and \\(\beta_{o}\\)*.
- **Consistency Polynomial:** *\\(\\{g^{\beta _{\ell} \ell _{i}(s) + \beta _{r} r _{i}(s) + \beta _{o} o _{i}(s)}\\} _{i \in [d]}\\)*
- **Proof Key:** Provided to the prover
  - \\(\\{g^{\ell_i(s)},g^{r_i(s)},g^{o_i(s)}\\}_{i\in[d]}\\)
  - \\(\\{g^{\alpha_{\ell} \ell_i(s)},g^{\alpha_{r} r_i(s)},g^{\alpha_{o} o_i(s)}\\}_{i\in[d]}\\)
  - \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Verification Key:**
  - \\(g^{t(s)}\\)
  - \\(g^{\alpha_{\ell}}\\), \\(g^{\alpha_{r}}\\), \\(g^{\alpha_{o}}\\)
  - *\\(g^{\beta_{\ell}}, g^{\beta_{r}}, g^{\beta_{o}}\\)*
- After distribution, the original \\(s\\), \\(\alpha_{\ell}\\), \\(\alpha_r\\), and \\(\alpha_o\\) values are securely destroyed.

**Protocol (Proving)**

- Execute the program and get the witness vector \\(w\\).
- Compute the linear-combinations of polynomials
  - \\(\ell(x) = \sum_{i=1}^{d} w_i \ell_{i}(x)\\)
  - \\(r(x) = \sum_{i=1}^{d} w_i r_{i}(x)\\)
  - \\(o(x) = \sum_{i=1}^{d} w_i o_{i}(x)\\)
- Compute \\(h(x) = \frac{\ell(x) r(x) - o(x)}{t(x)}\\)
- Evaluate each polynomial at \\(s\\):
  - \\(g^{\ell(s)} = \prod^{d}_{i=1} (g^{\ell_i(s)})^{w_i} \\)
  - \\(g^{r(s)} = \prod^{d}_{i=1} (g^{r_i(s)})^{w_i} \\)
  - \\(g^{o(s)} = \prod^{d}_{i=1} (g^{o_i(s)})^{w_i} \\)
- Evaluate each shifted polynomial at \\(s\\):
  - \\(g^{\alpha_{\ell} \ell(s)} = \prod^{d}_{i=1} (g^{\alpha _{\ell} \ell_i(s)})^{w_i} \\)
  - \\(g^{\alpha_{r} r(s)} = \prod^{d}_{i=1} (g^{\alpha _{r} r_i(s)})^{w_i} \\)
  - \\(g^{\alpha_{o} o(s)} = \prod^{d}_{i=1} (g^{\alpha _{o} o_i(s)})^{w_i} \\)
- *Evaluate each consistency polynomial at \\(s\\):*
  - *\\(g^{z(s)} = \prod^{d}_{i=1} g^{\beta _{\ell} \ell _{i}(s) + \beta _{r} r _{i}(s) + \beta _{o} o _{i}(s)}\\)*
- Calculate \\(g^{h(s)}\\) from \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- Proof: \\((\\) \\(g^{\ell(s)}, g^{r(s)}, g^{o(s)}, g^{\alpha_{\ell} \ell(s)}, g^{\alpha_{r} r(s)}, g^{\alpha_{o} o(s)}, g^{h(s)},\\) *\\(g^{z(s)}\\)* \\()\\)

**Protocol (Verification)**

- Parse proof as \\((g^{\ell}, g^r, g^o, g^{\ell'}, g^{r'}, g^{o'}, g^{h}, g^{z})\\)
- Check polynomial restrictions
  - \\(e(g^{\ell}, g^{\alpha_{\ell}}) = e(g^{\ell'}, g)\\)
  - \\(e(g^{r}, g^{\alpha_{r}}) = e(g^{r'}, g)\\)
  - \\(e(g^{o}, g^{\alpha_{o}}) = e(g^{o'}, g)\\)
- Check variable consistency
  - *\\(e(g^{\ell}, g^{\beta _{\ell}}) \cdot e(g^{r}, g^{\beta _{r}}) \cdot e(g^{o}, g^{\beta _{o}}) = e(g^{z}, g)\\)*
- Validity check
  - \\(e(g^{\ell}, g^{r}) = e(g^t,g^h) \cdot e(g^o, g)\\)

