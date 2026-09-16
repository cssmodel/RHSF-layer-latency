# RHSF-layer-latency

Data and Fortran code for Appendix C and Fig. 3 of

C. Shin and J. H. Chon, "Recursive Heaviside Dynamics and Nested-Sigmoid Gating: A Derivative-Cascade View of Layer Latency," *IEEE Access* (manuscript Access-2026-38113).

A full description is in `README.docx`.

## 1. What the code computes

The nest has 841 layers on the normalised grid tau_j = j/841, j = 1, ..., 841, so that tau_841 = 1 closes the nest. For each present layer k = 840, ..., 1 the present is substituted first, t = tau_k, and the state is evaluated:

```
G_841 = -t + tau_841
G_j   = -t + tau_j * H_{s_(j+1)}(G_(j+1)),   j = 840, ..., k
u_k   = 1 - H_{s_k}(G_k),   H_s(x) = 1 / (1 + exp(-2x/s))
```

The broad-gate sequence was drawn uniformly from [0.1, 11.1); the narrow-gate sequence is the same draw divided by 10, in [0.01, 1.11). The closing layer carries s_841 = 5 and 0.5 respectively. All arithmetic in `rhsf_from_file` is quadruple precision (kind=16).

## 2. Files

| File | Content |
|---|---|
| `real_high_s.txt` | Broad-gate input. 841 lines; column 1 = tau_j = j/841, column 2 = s_j; rows ordered j = 841, ..., 1. |
| `real_low_s.txt` | Narrow-gate input. Same layout; s_j is the broad-gate value divided by 10. |
| `tau.txt` | Column 1 = grid before normalisation (not used), column 2 = normalised tau_j. |
| `rhsf_from_file.f90` | Recomputes u_k from an input file. No random-number generator. |
| `real_high_s_u.asc`, `real_low_s_u.asc` | Reference output: column 1 = tau_k, column 2 = u_k, k = 840, ..., 1. |
| `high_s_unlimited.f90`, `low_s_unlimited.f90` | Original programs that drew s_j and wrote the input files. |

## 3. How to reproduce Fig. 3

```
gfortran -O2 rhsf_from_file.f90 -o rhsf_from_file
./rhsf_from_file real_high_s.txt real_high_s_u.asc
./rhsf_from_file real_low_s.txt  real_low_s_u.asc
```

Fig. 3(a), (b): column 2 against column 1 of `real_low_s.txt`, `real_high_s.txt`.
Fig. 3(c), (d): column 2 against column 1 of `real_low_s_u.asc`, `real_high_s_u.asc`.
The reference output was produced with gfortran 13.3.

## 4. Notes on reproducibility

- The input files, not the random seed, define the computation. The original programs set `seed = 20260730`, but the sequence `random_number` returns depends on the compiler and its version. Start from the `.txt` files.
- s_j was drawn in single precision and is stored to 8 significant digits; every later operation is in quadruple precision.
- An independent recomputation from the input files agrees with the reference output to within 1e-7.

## 5. License and contact

MIT License. Changsoo Shin, Department of Energy Resources Engineering, Seoul National University, cssmodel@snu.ac.kr.
