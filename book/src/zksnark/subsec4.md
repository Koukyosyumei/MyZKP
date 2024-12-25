# Bringing It All Together

Let's recap the previous sections. First, the relationship between the inputs and outputs of any program can be expressed as a rank-one constraint system (R1CS) as follows:

\\[
  (L \cdot v) \circ (R \cdot v) - O \cdot v = 0  
\\]

, where \\(v\\) is the concatenation of all inputs, outputs, and intermediate values. This allows us to transform the statement, "I know the input values \\(x\\) that make the program returns the output values \\(y\\)", into "I know \\(v\\), whose outputs components are \\(y\\), that satisfies the constraint system corresponding to the program". 

Then, instead of separately checking each constraint (which corresponds to a row in the R1CS matrix), we can convert this into a more efficient polynomial-equivalence test.

## First Protocol Idea

The simplest protocol, based on the previous chapter, is as follows:

**Protocol (Setup)**

- **Interpolated Polynomial:** Construct \\(\\{\ell_i, r_i, o_i\\}_{i\in[d]}\\) from \\(L\\), \\(R\\), and \\(O\\), respectively.
- **Target Polynomial:** \\(t(x) = (x-1)(x-2) \cdots (x-m)\\)
- **Secret Seed:** A trusted party generates the random value \\(s\\) and \\(\alpha\\).
- **Proof Key:** Provided to the prover
  - \\(\\{g^{\ell_i(s)},g^{r_i(s)},g^{o_i(s)}\\}_{i\in[d]}\\)
  - \\(\\{g^{\alpha \ell_i(s)},g^{\alpha r_i(s)},g^{\alpha o_i(s)}\\}_{i\in[d]}\\)
  - \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Verification Key:**
  - \\(g^{t(s)}, g^{\alpha}\\)
- After distribution, the original \\(s\\) and \\(\alpha\\) values are securely destroyed.

Both the proof key and the verification key are publicly available, enabling anyone to generate and verify proofs based on the target program.

**Protocol (Proving)**

- Run the program to obtain the assignment vector \\(v\\).
- Compute the linear-combinations of polynomials
  - \\(\ell(x) = \sum_{i=1}^{d} v_i \ell_{i}(x),\quad r(x) = \sum_{i=1}^{d} v_i r_{i}(x),\quad o(x) = \sum_{i=1}^{d} v_i o_{i}(x)\\)
- Compute the quotient polynomial:
  - \\(h(x) = \frac{\ell(x) r(x) - o(x)}{t(x)}\\)
- Evaluate each polynomial at \\(s\\).
  - \\(g^{\ell(s)} = \prod^{d}_{i=1} (g^{\ell_i(s)})^{v _i} ,\quad g^{r(s)} = \prod^{d} _{i=1} (g^{r_i(s)})^{v_i} ,\quad g^{o(s)} = \prod^{d} _{i=1} (g^{o _i(s)})^{v _i} \\)
- Evaluate the shifted polynomials at \\(s\\).
  - \\(g^{\alpha \ell(s)} = \prod^{d} _{i=1} (g^{\alpha \ell _i(s)})^{v _i} ,\quad g^{\alpha r(s)} = \prod^{d} _{i=1} (g^{\alpha r _i(s)})^{v _i} ,\quad g^{\alpha o(s)} = \prod^{d} _{i=1} (g^{\alpha o_i(s)})^{v _i} \\)
- Compute \\(g^{h(s)}\\) from \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Proof**: 
  - \\((g^{\ell(s)}, g^{r(s)}, g^{o(s)}, g^{\alpha \ell(s)}, g^{\alpha r(s)}, g^{\alpha o(s)}, g^{h(s)})\\)

**Protocol (Verification)**

- Parse the proof as \\((g^{\ell}, g^r, g^o, g^{\ell'}, g^{r'}, g^{o'}, g^{h})\\)
- Check the polynomial restrictions
  - \\(e(g^{\ell}, g^{\alpha}) = e(g^{\ell'}, g),\quad e(g^{r}, g^{\alpha}) = e(g^{r'}, g),\quad e(g^{o}, g^{\alpha}) = e(g^{o'}, g)\\)
- Verify validity of the proof
  - \\(e(g^{\ell}, g^{r}) = e(g^t,g^h) \cdot e(g^o, g)\\)

**Vulnerability**

A critical issue with this protocol is that the checks may pass even if \\(g^{\ell}\\) is computed not from \\(\\{g^{\ell_i(s)}\\}_{i \in [d]}\\) but from \\(\\{g^{r_i(s)}\\} _{i \in [d]}\\), \\(\\{g^{o_i(s)}\\} _{i \in [d]}\\), or their combinations. The same issue applies to \\(g^{r}\\) and \\(g^{o}\\). 

For example, if the prover sends \\((g^{r(s)}, g^{\ell(s)}, g^{o(s)}, g^{\alpha r(s)}, g^{\alpha \ell(s)}, g^{\alpha o(s)}, g^{h(s)})\\) as the proof, all the verification checks still pass, even although the proved statement differs from the original one.

## Second Protocol: Non-Interchangibility

To address the interchangeability issue, the next protocol uses distinct the different \\(\alpha\\)-shift for \\(\ell\\), \\(r\\), and \\(o\\).

**Protocol (Setup)**

- **Interpolated Polynomial:** Construct \\(\\{\ell_i, r_i, o_i\\}_{i\in[d]}\\) from \\(L\\), \\(R\\), and \\(O\\), respectively.
- **Target Polynomial:** \\(t(x) = (x-1)(x-2) \cdots (x-m)\\)
- **Secret Seed:** A trusted party generates the random value \\(s\\), *\\(\alpha_{\ell}\\), \\(\alpha_r\\), and \\(\alpha_o\\)*.
- **Proof Key (for the prover):**
  - \\(\\{g^{\ell_i(s)},g^{r_i(s)},g^{o_i(s)}\\}_{i\in[d]}\\)
  - *\\(\\{g^{\alpha_{\ell} \ell_i(s)},g^{\alpha_{r} r_i(s)},g^{\alpha_{o} o_i(s)}\\}_{i\in[d]}\\)*
  - \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Verification Key (public):**
  - \\(g^{t(s)}\\) *\\(, g^{\alpha_{\ell}},g^{\alpha_{r}},g^{\alpha_{o}}\\)*
- After distribution, the original \\(s\\), \\(\alpha_{\ell}\\), \\(\alpha_r\\), and \\(\alpha_o\\) values are securely destroyed.

**Protocol (Proving)**

- Execute the program and get the assignment \\(v\\).
- Compute the linear-combinations of polynomials
  - \\(\ell(x) = \sum_{i=1}^{d} v_i \ell_{i}(x),\quad r(x) = \sum_{i=1}^{d} v_i r_{i}(x),\quad o(x) = \sum_{i=1}^{d} v_i o_{i}(x)\\)
- Compute the quotient polynomial:
  - \\(h(x) = \frac{\ell(x) r(x) - o(x)}{t(x)}\\)
- Evaluate each polynomial at \\(s\\).
  - \\(g^{\ell(s)} = \prod^{d}_{i=1} (g^{\ell_i(s)})^{v _i} ,\quad g^{r(s)} = \prod^{d} _{i=1} (g^{r_i(s)})^{v_i} ,\quad g^{o(s)} = \prod^{d} _{i=1} (g^{o _i(s)})^{v _i} \\)
- Evaluate each shifted polynomial at \\(s\\).
  - *\\(g^{\alpha_{\ell} \ell(s)} = \prod^{d}_{i=1} (g^{\alpha _{\ell} \ell_i(s)})^{v_i} \\)*, *\\(g^{\alpha_{r} r(s)} = \prod^{d}_{i=1} (g^{\alpha _{r} r_i(s)})^{v_i} \\)*, *\\(g^{\alpha_{o} o(s)} = \prod^{d}_{i=1} (g^{\alpha _{o} o_i(s)})^{v_i} \\)*
- Calculate \\(g^{h(s)}\\) from \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Proof**: 
  - \\((g^{\ell(s)}, g^{r(s)}, g^{o(s)},\\) *\\(g^{\alpha_{\ell} \ell(s)}, g^{\alpha_{r} r(s)}, g^{\alpha_{o} o(s)},\\)* \\(g^{h(s)})\\)

**Protocol (Verification)**

- Parse proof as \\((g^{\ell}, g^r, g^o, g^{\ell'}, g^{r'}, g^{o'}, g^{h})\\)
- Check polynomial restrictions
  - *\\(e(g^{\ell}, g^{\alpha_{\ell}}) = e(g^{\ell'}, g)\\)*, \\(\quad \\) *\\(e(g^{r}, g^{\alpha_{r}}) = e(g^{r'}, g)\\)*, \\(\quad \\) *\\(e(g^{o}, g^{\alpha_{o}}) = e(g^{o'}, g)\\)*
- Validity check
  - \\(e(g^{\ell}, g^{r}) = e(g^t,g^h) \cdot e(g^o, g)\\)

**Vulnerability**

This protocol resolves interchangeability but does not enforce consistency acros \\(\ell\\), \\(r\\), and \\(o\\). Variables \\(v_i\\) can still take different values in each polynoimal because verification checks are performed separately.

## Thrid Protocol: Variable Consistency

To achive the variable consistency, we employ a checksum mechanism. Specifically, we first draw a new random value \\(\beta\\) and define the checksum of \\(v_i\\) as \\(g^{\beta(\ell_{i}(s) + r_i(s) + o_i(s))}\\). Let \\(v _{\ell, i}\\), \\(v _{r,i}\\), \\(v _{o,i}\\), and \\(v _{\beta,i}\\) denote the \\(i\\)-th value of the assignment vectors for \\(\ell\\), \\(r\\), \\(o\\) and the checksum, respectively. If all of them are the same, the following equation holds:

\\[
e(g^{v _{\ell, i} \ell_i(s)} g^{v _{r, i} r_i(s)} g^{v _{o, i} o_i(s)}, g^{\beta}) = e(g^{v _{\beta, i} \beta(\ell _{i}(s) + r _{i}(s) + o _{i}(s))}, g)
\\]

Unfortunately, this condition is not strictly equivalent. For example, consider the case where \\(\ell_i(x) = r_i(x)\\). In this scenario, we have:

\begin{align*}
  &\beta(v _{\ell, i} \ell_i(s) + v _{r, i} r_i(s) + v _{o, i} o_i(s)) = v _{\beta, i} \beta (\ell _{i}(s) + r _{i}(s) + o _{i}(s)) \\\\
  \iff &\beta(v _{\ell, i} \ell_i(s) + v _{r, i} \ell_i(s) + v _{o, i} o_i(s)) = v _{\beta, i} \beta (2\ell _{i}(s) + o _{i}(s))
\end{align*}

This equation holds for arbitrary \\(v _{r,i}\\) and \\(v _{o,i}\\) if we set \\(v _{\beta, i} = v _{o, i}\\) and \\(v _{\ell, i} = 2 v _{o, i} - v _{r, i}\\).

To address this issue, we use distinct different \\(\beta\\) values for \\(\ell\\), \\(r\\) and \\(o\\). The consistency check then verifies the following equation:

\\[
e(g^{v _{\ell, i} \ell _{i}(s)}, g^{\beta _{\ell}}) \cdot e(g^{v _{r, i} r _{i}(s)}, g^{\beta _{r}}) \cdot e(g^{v _{o, i} o _{i}(s)}, g^{\beta _{o}}) = e(g^{v _{\beta, i}(\beta _{\ell} \ell _{i}(s) + \beta _{r} r _{i}(s) + \beta _{o} o _{i}(s))}, g)   
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
  - \\(g^{t(s)}, g^{\alpha_{\ell}},g^{\alpha_{r}},g^{\alpha_{o}}\\) *\\(, g^{\beta_{\ell}}, g^{\beta_{r}}, g^{\beta_{o}}\\)*
- After distribution, the original \\(s\\), \\(\alpha_{\ell}\\), \\(\alpha_r\\), and \\(\alpha_o\\) values are securely destroyed.

**Protocol (Proving)**

- Execute the program and get the assignment vector \\(v\\).
- Compute the linear-combinations of polynomials
  - \\(\ell(x) = \sum_{i=1}^{d} v_i \ell_{i}(x),\quad r(x) = \sum_{i=1}^{d} v_i r_{i}(x),\quad o(x) = \sum_{i=1}^{d} v_i o_{i}(x)\\)
- Compute the quotient polynomial:
  - \\(h(x) = \frac{\ell(x) r(x) - o(x)}{t(x)}\\)
- Evaluate each polynomial at \\(s\\):
  - \\(g^{\ell(s)} = \prod^{d}_{i=1} (g^{\ell_i(s)})^{v _i} ,\quad g^{r(s)} = \prod^{d} _{i=1} (g^{r_i(s)})^{v_i} ,\quad g^{o(s)} = \prod^{d} _{i=1} (g^{o _i(s)})^{v _i} \\)
- Evaluate each shifted polynomial at \\(s\\):
  - \\(g^{\alpha _{\ell} \ell(s)} = \prod^{d} _{i=1} (g^{\alpha _{\ell} \ell _i(s)})^{v _i} ,g^{\alpha _{r} r(s)} = \prod^{d} _{i=1} (g^{\alpha _{r} r _i(s)})^{v _i} ,g^{\alpha _{o} o(s)} = \prod^{d} _{i=1} (g^{\alpha _{o} o _i(s)})^{v _i} \\)
- *Evaluate each consistency polynomial at \\(s\\):*
  - *\\(g^{z(s)} = \prod^{d}_{i=1} (g^{\beta _{\ell} \ell _{i}(s) + \beta _{r} r _{i}(s) + \beta _{o} o _{i}(s)})^{v _{i}}\\)*
- Calculate \\(g^{h(s)}\\) from \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Proof**: 
  - \\((\\) \\(g^{\ell(s)}, g^{r(s)}, g^{o(s)}, g^{\alpha_{\ell} \ell(s)}, g^{\alpha_{r} r(s)}, g^{\alpha_{o} o(s)}, g^{h(s)},\\) *\\(g^{z(s)}\\)* \\()\\)

**Protocol (Verification)**

- Parse proof as \\((g^{\ell}, g^r, g^o, g^{\ell'}, g^{r'}, g^{o'}, g^{h}, g^{z})\\)
- Check polynomial restrictions
  - \\(\quad e(g^{\ell}, g^{\alpha_{\ell}}) = e(g^{\ell'}, g),\quad e(g^{r}, g^{\alpha_{r}}) = e(g^{r'}, g),\quad e(g^{o}, g^{\alpha_{o}}) = e(g^{o'}, g)\\)
- Check variable consistency
  - *\\(e(g^{\ell}, g^{\beta _{\ell}}) \cdot e(g^{r}, g^{\beta _{r}}) \cdot e(g^{o}, g^{\beta _{o}}) = e(g^{z}, g)\\)*
- Validity check
  - \\(e(g^{\ell}, g^{r}) = e(g^t,g^h) \cdot e(g^o, g)\\)

**Vulnerability**

Despite these checks, the protocol is vulnerable to **malleability**. Specifically, a malicious prover can exploit the polynomial restriction check to introduce arbitrary constants, altering the proof witout detection.

Recall that the verifier validates whether the submitted \\(g^{\ell}\\) is actually calculated by \\(\\{g^{\ell _i(s)}\\} _{i \in [d]}\\) by checking \\(e(g^{\ell}, g^{\alpha _{\ell}}) = e(g^{\ell'}, g)\\). However, this process is not sound. Recall that the verification key is publicaly available, and the prover knows both of \\(g^{\alpha _{\ell}}\\) and \\(g^{\beta _{\ell}}\\). Here, suppose the prover submits \\(g^{\ell} g^{c}\\) and \\(g^{\ell'} (g^{\alpha _{\ell}})^{c}\\) insteads of \\(g^{\ell}\\) and \\(g^{\ell'}\\), where \\(c\\) is a constatn value. Then, the polynomial restriction check still passes:

\\[
e(g^{\ell} g^{c}, g^{\alpha _{\ell}}) = e(g^{\alpha _{\ell} \ell + \alpha _{\ell} c}, g) = e(g^{\ell'}g^{\alpha _{\ell}c}, g) = e(g^{\ell'} (g^{\alpha _{\ell}})^{c}, g)
\\]

In addition, if the prover submits \\(g^{z} (g^{\beta _{\ell}})^{c}\\) as the checksum, it also passes the polynomial checksum verification:

\\[
e(g^{\ell} g^{c}, g^{\beta _{\ell}}) \cdot e(g^{r}, g^{\beta _{r}}) \cdot e(g^{o}, g^{\beta _{o}}) = e(g^{z} (g^{\beta _{\ell}})^{c}, g)
\\]

This phenomenon also can occur for \\(r\\) and \\(o\\).

## Forth Protocol: Non-Malleability

One way to surrogate the above malleability is hiding \\(g^{\beta _{\ell}}\\), \\(g^{\beta _{r}}\\), and \\(g^{\beta _{o}}\\) by powering them with a new random value \\(\eta\\).


**Protocol (Setup)**

- **Interpolated Polynomial:** Construct \\(\\{\ell_i, r_i, o_i\\}_{i\in[d]}\\) from \\(L\\), \\(R\\), and \\(O\\), respectively.
- **Target Polynomial:** \\(t(x) = (x-1)(x-2) \cdots (x-m)\\)
- **Secret Seed:** A trusted party generates the random value \\(s\\), \\(\alpha_{\ell}\\), \\(\alpha_r\\), \\(\alpha_o\\), \\(\beta_{\ell}\\), \\(\beta_{r}\\), \\(\beta_{o}\\), and *\\(\eta\\)*.
- **Consistency Polynomial:** \\(\\{g^{\beta _{\ell} \ell _{i}(s) + \beta _{r} r _{i}(s) + \beta _{o} o _{i}(s)}\\} _{i \in [d]}\\)
- **Proof Key:** Provided to the prover
  - \\(\\{g^{\ell_i(s)},g^{r_i(s)},g^{o_i(s)}\\}_{i\in[d]}\\)
  - \\(\\{g^{\alpha_{\ell} \ell_i(s)},g^{\alpha_{r} r_i(s)},g^{\alpha_{o} o_i(s)}\\}_{i\in[d]}\\)
  - \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Verification Key:**
  - \\(g^{t(s)}, g^{\alpha_{\ell}},g^{\alpha_{r}},g^{\alpha_{o}}\\) *\\(, g^{\beta_{\ell} \eta}, g^{\beta_{r} \eta}, g^{\beta_{o} \eta}\\)*
- After distribution, the original \\(s,\alpha_{\ell},\alpha_r\\), and \\(\alpha_o\\) values are securely destroyed.

**Protocol (Proving)**

- Execute the program and get the assignment vector \\(v\\).
- Compute the linear-combinations of polynomials
  - \\(\ell(x) = \sum_{i=1}^{d} v_i \ell_{i}(x),\quad r(x) = \sum_{i=1}^{d} v_i r_{i}(x),\quad o(x) = \sum_{i=1}^{d} v_i o_{i}(x)\\)
- Compute the quotient polynomial:
  - \\(h(x) = \frac{\ell(x) r(x) - o(x)}{t(x)}\\)
- Evaluate each polynomial at \\(s\\):
  - \\(g^{\ell(s)} = \prod^{d}_{i=1} (g^{\ell_i(s)})^{v _i} ,\quad g^{r(s)} = \prod^{d} _{i=1} (g^{r_i(s)})^{v_i} ,\quad g^{o(s)} = \prod^{d} _{i=1} (g^{o _i(s)})^{v _i} \\)
- Evaluate each shifted polynomial at \\(s\\):
  - \\(g^{\alpha _{\ell} \ell(s)} = \prod^{d} _{i=1} (g^{\alpha _{\ell} \ell _i(s)})^{v _i} ,g^{\alpha _{r} r(s)} = \prod^{d} _{i=1} (g^{\alpha _{r} r _i(s)})^{v _i} ,g^{\alpha _{o} o(s)} = \prod^{d} _{i=1} (g^{\alpha _{o} o _i(s)})^{v _i} \\)
- Evaluate each consistency polynomial at \\(s\\):
  - \\(g^{z(s)} = \prod^{d}_{i=1} (g^{\beta _{\ell} \ell _{i}(s) + \beta _{r} r _{i}(s) + \beta _{o} o _{i}(s)})^{v _{i}}\\)
- Calculate \\(g^{h(s)}\\) from \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Proof**: 
  - \\((\\) \\(g^{\ell(s)}, g^{r(s)}, g^{o(s)}, g^{\alpha_{\ell} \ell(s)}, g^{\alpha_{r} r(s)}, g^{\alpha_{o} o(s)}, g^{h(s)},\\) \\(g^{z(s)}\\) \\()\\)

**Protocol (Verification)**

- Parse proof as \\((g^{\ell}, g^r, g^o, g^{\ell'}, g^{r'}, g^{o'}, g^{h}, g^{z})\\)
- Check polynomial restrictions
  - \\(e(g^{\ell}, g^{\alpha_{\ell}}) = e(g^{\ell'}, g),\quad e(g^{r}, g^{\alpha_{r}}) = e(g^{r'}, g),\quad e(g^{o}, g^{\alpha_{o}}) = e(g^{o'}, g)\\)
- Check variable consistency
  - *\\(e(g^{\ell}, g^{\beta _{\ell} \eta}) \cdot e(g^{r}, g^{\beta _{r} \eta}) \cdot e(g^{o}, g^{\beta _{o} \eta}) = e(g^{z}, g^{\eta})\\)*
- Validity check
  - \\(e(g^{\ell}, g^{r}) = e(g^t,g^h) \cdot e(g^o, g)\\)

## Fifth Protocol: Pinocchio

The above protocol introduces four expensive pairing operations. To make it faster, Pinocchio protocl randomizes the generators.

**Protocol (Setup)**

- **Interpolated Polynomial:** Construct \\(\\{\ell_i, r_i, o_i\\}_{i\in[d]}\\) from \\(L\\), \\(R\\), and \\(O\\), respectively.
- **Target Polynomial:** \\(t(x) = (x-1)(x-2) \cdots (x-m)\\)
- **Secret Seed:** A trusted party generates the random value \\(s\\), \\(\alpha_{\ell}\\), \\(\alpha_r\\), \\(\alpha_o\\), \\(\beta_{\ell}\\), \\(\beta_{r}\\), \\(\beta_{o}\\), \\(\eta\\), *\\(\rho _{\ell}\\), and \\(\rho _{r}\\), and set \\(\rho _{o} = \rho _{\ell} \rho _{r}\\)*.
- **Randized Generators:** *\\(g _{\ell} = g^{\rho _{\ell}}\\), \\(g _{r} = g^{\rho _{r}}\\), and \\(g _{o} = g^{\rho _{o}}\\)*
- **Consistency Polynomial:** \\(\\{g^{\beta _{\ell} \ell _{i}(s) + \beta _{r} r _{i}(s) + \beta _{o} o _{i}(s)}\\} _{i \in [d]}\\)
- **Proof Key:** Provided to the prover
  - *\\(\\{g_{\ell}^{\ell_i(s)},g_{r}^{r_i(s)},g_{o}^{o_i(s)}\\}_{i\in[d]}\\)*
  - *\\(\\{g_{\ell}^{\alpha_{\ell} \ell_i(s)},g_{r}^{\alpha_{r} r_i(s)},g_{o}^{\alpha_{o} o_i(s)}\\}_{i\in[d]}\\)*
  - \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Verification Key:**
  - *\\(g _{o}^{t(s)}\\)* \\(, g^{\alpha _{\ell}}, g^{\alpha _{r}}, g^{\alpha _{o}}, g^{\beta _{\ell} \eta}, g^{\beta _{r} \eta}, g^{\beta _{o} \eta}\\)
- After distribution, the original \\(s\\), \\(\alpha_{\ell}\\), \\(\alpha_r\\), and \\(\alpha_o\\) values are securely destroyed.

**Protocol (Proving)**

- Execute the program and get the assignment vector \\(v\\).
- Compute the linear-combinations of polynomials
  - \\(\ell(x) = \sum_{i=1}^{d} v_i \ell_{i}(x),\quad r(x) = \sum_{i=1}^{d} v_i r_{i}(x),\quad o(x) = \sum_{i=1}^{d} v_i o_{i}(x)\\)
- Compute the quotient polynomial:
  - \\(h(x) = \frac{\ell(x) r(x) - o(x)}{t(x)}\\)
- Evaluate each polynomial at \\(s\\):
  - \\(g^{\ell(s)} = \prod^{d}_{i=1} (g^{\ell_i(s)})^{v _i} ,\quad g^{r(s)} = \prod^{d} _{i=1} (g^{r_i(s)})^{v_i} ,\quad g^{o(s)} = \prod^{d} _{i=1} (g^{o _i(s)})^{v _i} \\)
- Evaluate each shifted polynomial at \\(s\\):
  - \\(g^{\alpha _{\ell} \ell(s)} = \prod^{d} _{i=1} (g^{\alpha _{\ell} \ell _i(s)})^{v _i} ,g^{\alpha _{r} r(s)} = \prod^{d} _{i=1} (g^{\alpha _{r} r _i(s)})^{v _i} ,g^{\alpha _{o} o(s)} = \prod^{d} _{i=1} (g^{\alpha _{o} o _i(s)})^{v _i} \\)
- Evaluate each consistency polynomial at \\(s\\):
  - *\\(g^{z(s)} = \prod^{d}_{i=1} (g _{\ell}^{\beta _{\ell} \ell _{i}(s)} \cdot g _{r}^{\beta _{r} r _{i}(s)} \cdot g _{o}^{\beta _{o} o _{i}(s)})^{v _{i}}\\)*
- Calculate \\(g^{h(s)}\\) from \\(\\{g^{(s^j)}\\}_{j\in[m]}\\)
- **Proof**: 
  - \\((\\) \\(g^{\ell(s)}, g^{r(s)}, g^{o(s)}, g^{\alpha_{\ell} \ell(s)}, g^{\alpha_{r} r(s)}, g^{\alpha_{o} o(s)}, g^{h(s)},\\) \\(g^{z(s)}\\) \\()\\)

**Protocol (Verification)**

- Parse proof as \\((g^{\ell}, g^r, g^o, g^{\ell'}, g^{r'}, g^{o'}, g^{h}, g^{z})\\)
- Check polynomial restrictions
  - *\\(e(g _{\ell}^{\ell}, g^{\alpha _{\ell}}) = e(g _{\ell}^{\ell'}, g)\\)*, *\\(e(g _{r}^{r}, g^{\alpha _{r}}) = e(g _{r}^{r'}, g)\\)*, *\\(e(g _{o}^{o}, g^{\alpha _{o}}) = e(g _{o}^{o'}, g)\\)*
- Check variable consistency
  - *\\(e(g _{\ell}^{\ell} \cdot g _{r}^{r} \cdot g _{o}^{o}) = e(g^{z}, g^{\eta})\\)*
- Validity check
  - \\(e(g^{\ell}, g^{r}) = e(g^t,g^h) \cdot e(g^o, g)\\)

This protocol reduces two pairing operations.