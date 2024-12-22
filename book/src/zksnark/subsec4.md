# Bringing It All Together

We first recap the previous sections. First the relationship between inputs and outputs of any prgram can be represented as a rank-one constraint system (R1CS) as follows:

\\[
  (L \cdot a) \circ (R \cdot a) - O \cdot a = 0  
\\]

, where \\(a\\) is the concatenation of all inputs, outputs, and intermediate values (witness). Then, the statement, "I know the input values \\(x\\) that make the program returns the output values \\(y\\)", can be converted into "I know \\(a\\), whose outputs part is \\(y\\), that satisfies the constraint system corresponding to the program". 

Although separately checking each constraint, corresponding to each row in the above R1CS, is not efficient, we can convert those vector-equivalence test into the polynomial-equivalence test.

## 