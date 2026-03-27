# README

This folder contains the Python scripts used to compute solutions of Eq. (33), verify the correction orders and nonlinearities of the examples presented in the paper, and validate the derivations of several coefficients used in the theoretical analysis.

## File descriptions

### `Search solutions for m=2.py`
This script is used to solve Eq. (33) in the case \(m=2\) and \(t=n-m-1\).  
The program outputs the solutions for \(n=3,4,5,6\). The corresponding results are presented in **Example III.1** of the paper.

### `Search solutions for m=3.py`
This script is used to solve Eq. (33) in the case \(m=3\) and \(t=n-m-1\).  
The program outputs the solutions for \(n\in[4,12]\). The corresponding results are used in **Example V.1** of the paper.

### `Symmetric correctors.py`
This script is used to solve Eq. (33) in the case \(m=1\), under the assumption that for every \(0\le j\le n-t-1\),
\[
Y_{i,j}=\binom{n}{j}\quad \text{or}\quad 0.
\]
The corresponding results are used in **Remark III.1** of the paper. These solutions give the weight distributions associated with symmetric correctors.

### `Verify Example III.1.py`
This script is used to verify the correction order and nonlinearity of the function \(F(x)\) given in **Example III.1**.

### `Verify Example III.3.py`
This script is used to verify the correction orders and nonlinearities of the functions \(f_M\) and \(f_{M'}\) given in **Example III.3**.

### `Verify Example V.1.py`
This script is used to verify the correction orders and nonlinearities of the functions \(F(x)\) and \(H(F(x))\) given in **Example V.1**.

### `Verify Example V.3.py`
This script is used to verify the correction order and nonlinearity of the function \(F(x)\) given in **Example V.3**.

### `Verify St and St1.py`
This script is used to verify that the formulas for the parameters \(S_t\) and \(S_{t-1}\) in **Lemma IV.2** hold for all \(n\le 20\).

### `Verify Tjt and Tjt1.py`
This script is used to verify that the formulas for the parameters \(T_{j,t}\) and \(T_{j,t-1}\) in **Lemma IV.1** hold for all \(n\le 20\).

## Notes
- The scripts are intended to support the computational results appearing in the examples and remarks of the paper.
- The notation used in the scripts is consistent with that in the manuscript.
