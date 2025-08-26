# Gemini

While many polynomial commitment schemes work for standard single-variable (univariate) polynomials, things get trickier when we deal with multilinear polynomials, which are polynomials with multiple variables where each variable has a maximum degree of one. The Gemini protocol offers an elegant way to handle these.

## Multilinear Polynomial

First, let's get a handle on the main object of interest. A multilinear polynomial with \\(m\\) variables (\\(f(X_0, X_1, \dots, X_{m- 1})\\)) can be written as a sum over all possible combinations of its variables. A general \\(m\\)-variate multilinear polynomial \\(f\\) is defined as:

\begin{align*}
    f(X_0, X_1, \dots, X_{m-1}) = \sum_{i = 0}^{2^m-1} c_i \prod_{j=0}^{m-1} X_j^{b_j(i)}
\end{align*}

Here, \\(c_i\\) is the \\(i\\)-th coefficient, and \\((b_{m-1}(i), \dots, b_0(i))\\) is the binary representation of the number \\(i\\). THis formula just means we sum up terms where each variable \\(X_{j}\\) is either present or absent.

For example, a three-variable (\\(m = 3\\)) multilinear polynomial looks like this:

\begin{align*}
    f(X_0, X_1, X_2) = c_{000} + c_{100} X_0 + c_{010} X_1 + c_{110} X_0 X_1 + c_{001} X_2 + c_{101} X_0 X_2 + c_{011} X_1 X_2 + c_{111} X_0 X_1 X_2
\end{align*}

The subscript of each coefficient \\(c\\) tells you which variables are in that term. For instance, \\(c_{110}\\) is the coefficient of the \\(X_0 X_1\\) term (reading the bits from right to left as \\(i_0\\),\\(i_1\\),\\(i_2\\)).

## Problem

Our goal is to create a commitment scheme for these polynomials. Specifically, a Prover wants to commit to a multilinear polynomial \\(f\\) and then convince a Verifier that for a specific, publicly known point \\(\rho = (\rho_0, \rho_1, \dots, \rho_{m-1})\\) and a public value \\(u\\), the following statement is true:

\begin{align*}
    f(\rho) = u
\end{align*}

The Prover must prove this without revealing the entire polynomial \\(f\\).

## Key idea 

Here's the core idea of Gemini: **transform the multivariable problem into a series of single-variable problems**, which are much easier to handle.

We first take the \\(2^{m}\\) coefficients of our multilinear polynomial \\(f\\) and use them to define a new **univariate** polynomial \\(g\\) of degree \\(2^{m} - 1\\). 

\begin{align*}
    g(X) = \sum_{i=0}^{2^{m} - 1} c_{i} X^{i}
\end{align*}

Using our three-variable example, the corresponding univariate polynomial \\(g(X)\\) would be:

\begin{align*}
g(X) = c_{000} + c_{100} X + c_{010} X^2 + c_{110} X^3 + c_{001} X^4 + c_{101} X^5 + c_{011} X^6 + c_{111} X^7
\end{align*}

### Split and Fold

Then, we are going to recursively fold this univariate polynomialan univariate polynomial with the **folding** trick, which repeatedly halves the degree by separating even and odd exponents:

For any univariate polynomial \\(g(X)\\), we can write:

\begin{align*}
    g(X) = g_E(X) + X \cdot g_O(X)
\end{align*}, where

\begin{align*}
    g_E(X) = \frac{g(X) + g(-X)}{2}, \quad \quad g_O(X) = \frac{g(X) - g(-X)}{2X}
\end{align*}

Here, \\(g_E\\) contains the even exponents of \\(g\\) and \\(g_O\\) contains the shifted odd exponents.

\begin{align*}
    g_E(X) = \sum_{i=0}^{\frac{m+1}{2} - 1} c_{2i} X^{2i}, \quad \quad g_O(X) = \sum_{i=0}^{\frac{m+1}{2} - 1} c_{2i + 1} X^{2i}
\end{align*}

Now treat \\(Y := X^2\\). Because the even and odd parts only use powers \\(X^{2k}\\), both \\(g_E\\) and \\(g_O\\) can be seen as polynomials in \\(Y\\) with half the degree.

We then define a sequence of folded polynomials as follows:

\begin{align*}
    g^{(0)}(X) &:= g(X) \\\\
    g^{(1)}(X) &:= g_E^{(0)}(X) + \rho_0 g_O^{(0)}(X) \\\\
    g^{(2)}(X^2) &:= g_E^{(1)}(X^2) + \rho_1 g_O^{(1)}(X^2) \\\\
        &\quad\vdots \\\\
    g^{(m)}(X^{2^{m-1}}) &:= g_E^{(m-1)}(X^{2^{m-1}}) + \rho_{m-1} g_O^{(m-1)}(X^{2^{m-1}}) = (\text{a constant}) \\\\
\end{align*}

Intuitively, each folding step substitutes the known scalar \\(\rho_{i-1}\\) for one of the original variables, and reduces the number of remaining variables by one. After \\(m\\) folds, we obtain a constant equal to \\(f(\rho)\\), i.e., \\(g^{(m)}(\cdot) = u\\).

For example, the above three-variable example yeilds the following sequence, where \\(Y_i\\) denotes \\(X^{2^{i}}\\):

\begin{align*}
g^{(0)}(Y_0) &= c_{000} + c_{100} Y_0 + c_{010} Y_0^2 + c_{110} Y_0^3 + c_{001} Y_0^4 + c_{101} Y_0^5 + c_{011} Y_0^6 + c_{111} Y_0^7 \\\\
             &= c_{000} + c_{010} Y_0^2 + c_{001} Y_0^4 + c_{011} Y_0^6 + Y_0 (c_{100} + c_{110} Y_0^2 + c_{101} Y_0^4 + c_{111} Y_0^6) \\\\
g^{(1)}(Y_1) &= g_{E}^{(0)}(Y_1) + \rho_0 g_{O}^{(0)}(Y_1) \\\\ 
             &= c_{000} + c_{010} Y_1 + c_{001} Y_1^2 + c_{011} Y_1^3 + \rho_0 (c_{100} + c_{110} Y_1 + c_{101} Y_1^2 + c_{111} Y_1^3) \\\\
             &= (c_{000} + \rho_0  c_{100}) + (c_{010} + \rho_0  c_{110}) Y_1 + (c_{001} + \rho_0  c_{101}) Y_1^2 + (c_{011} + \rho_0  c_{111}) Y_1^3 \\\\
             &= (c_{000} + \rho_0  c_{100}) + (c_{001} + \rho_0  c_{101}) Y_1^2 + Y_1 ((c_{010} + \rho_0  c_{110}) + (c_{011} + \rho_0  c_{111}) Y_1^2) \\\\
g^{(2)}(Y_2) &= g_{E}^{(1)}(Y_2) + \rho_1 g_{O}^{(1)}(Y_2) \\\\ 
             &= (c_{000} + \rho_0  c_{100}) + (c_{001} + \rho_0  c_{101}) Y_2 + \rho_1 ((c_{010} + \rho_0  c_{110}) + (c_{011} + \rho_0  c_{111}) Y_2) \\\\
             &= (c_{000} + \rho_0  c_{100} + \rho_1 c_{010} + \rho_0 \rho_1 c_{110}) + Y_2 (c_{001} + \rho_0 c_{101} + \rho_1 c_{011} + \rho_0 \rho_1 c_{111}) \\\\
g^{(3)}(Y_3) &= g_{E}^{(1)}(Y_3) + \rho_2 g_{O}^{(1)}(Y_3) \\\\
             &= (c_{000} + \rho_0  c_{100} + \rho_1 c_{010} + \rho_0 \rho_1 c_{110}) + \rho_2 (c_{001} + \rho_0 c_{101} + \rho_1 c_{011} + \rho_0 \rho_1 c_{111}) \\\\
             &= c_{000} + \rho_0  c_{100} + \rho_1 c_{010} + \rho_0 \rho_1 c_{110} + \rho_2 c_{001} + \rho_0 \rho_2 c_{101} + \rho_1 \rho_2 c_{011} + \rho_0 \rho_1 \rho_2 c_{111} \\\\
             &= u
\end{align*}

### Relationship used for verification

The recurrence used between consecutive folded polynomials is:

\begin{align*}
    g^{(i)}(X^2) &= g_{E}^{(i-1)}(X) + \rho_{i-1} g_{O}^{(i-1)}(X) \\\\ 
               &=\frac{g^{(i-1)}(X) + g^{(i-1)}(-X)}{2} + \rho_{i-1} \frac{g^{(i-1)}(X) - g^{(i-1)}(-X)}{2X}
\end{align*}



## Proving

Then, to prove the original claim saying that "I know a multilinear polynomial \\(f\\) whose evaluation on \\(\vec{\rho}\\) is \\(u\\)", where \\(\vec{\rho}\\) and \\(u\\) are public, the verifier first draws a random value \\(\beta\\) asks the prover to open each \\(g^{i}\\) on \\(\beta\\), \\(-\beta\\) and \\(\beta^2\\). Then, the verifier checks the above relation between \\(g^{(i-1)}\\) and \\(g^{(i)}\\), and whether \\(g^{(m)}(\beta) = u\\). This verification of each polynomial can be implemented with KZG commitmentm, and evaluation on \\(\beta\\), \\(-\beta\\) and \\(\beta^2\\) can be done with a batch proof of KZG commitmment.

For example, suppose the following polynomial \\(f(X_0, X_1, X_2) = 1 + 2 X_0 + 3 X_1 + 4 X_0 X_1 + 5 X_2 + 6 X_0 X_2 + 7 X_1 X_2 + 8 X_0 X_1 X_2\\), and the following evaluatin point \\(\rho_0 = 1, \rho_1 = 2, \rho_2 = 3\\). \\(f(\rho_0, \rho_1, \rho_2) = 140\\). 

\begin{align*}
    g^{(0)}(Y_0) &= 1 + 3 Y_0^2 + 5 Y_0^4 + 7 Y_0^6 + Y_0 (2 + 4 Y_0^2 + 6 Y_0^4 + 8 Y_0^6)  \\\\
    g^{(1)}(Y_1) &= 1 + 3 Y_1 + 5 Y_1^2 + 7 Y_1^3 + (2 + 4 Y_1 + 6 Y_1^2 + 8 Y_1^3) \\\\
                 &= 3 + 7 Y_1 + 11 Y_1^2 + 15 Y_1^3\\\\
                 &= 3 + 11 Y_1^2 + Y_1 (7 + 15 Y_1^2) \\\\
    g^{(2)}(Y_2) &= 3 + 11 Y_2 + 2 (7 + 15 Y_2)\\\\
                 &= 17 + Y_2 (41) \\\\
    g^{(3)}(Y_3) &= 17 + 3 \cdot 41 = 140
\end{align*}

Let's use \\(\beta = 2\\).

\begin{align*}
    g^{(0)}(2) &= 1793, \quad g^{(0)}(-2) = -711, \quad g^{(1)}(4) = 1167, \quad \frac{g^{(0)}(2) + g^{(0)}(-2)}{2} + 1 \cdot \frac{g^{(0)}(2) - g^{(0)}(-2)}{2 * 2} = 1167 \\\\
    g^{(1)}(2) &= 181, \quad g^{(1)}(-2) = -87, \quad g^{(2)}(4) = 181, \quad \frac{g^{(1)}(2) + g^{(1)}(-2)}{2} + 2 \cdot \frac{g^{(1)}(2) - g^{(1)}(-2)}{2 * 2} = 181 \\\\
    g^{(2)}(2) &= 99, \quad g^{(2)}(-2) = -65, \quad g^{(3)}(4) = u = 140, \quad \frac{g^{(2)}(2) + g^{(2)}(-2)}{2} + 3 \cdot \frac{g^{(2)}(2) - g^{(2)}(-2)}{2 * 2} = 140
\end{align*}
 