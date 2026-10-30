# DARTS

This repository contains the revised specification, implementations, and test vectors for the DARTS lattice-based digital signature scheme.

## Update

Security analysis identified a bias in the distribution of `b` caused by a missing factor of 2 in Step 9 of Algorithm 7. This release changes the comparison from

`Comp_d(ω)` and `Comp_d(ω + (1 - b)cα)`

to

`Comp_d(ω)` and `Comp_d(ω + (1 - 2b)cα)`.

The correction has been applied to DARTS-128, DARTS-256, and DARTS-512, together with their test vectors. The average number of rejection-sampling attempts increases by less than one for each parameter set, while the DARTS-512 signature grows by one byte. The updated sizes are shown below.

| Parameter set | Public key (bytes) | Secret key (bytes) | Signature (bytes) |
|---|---:|---:|---:|
| DARTS-128 | 1,120 | 1,536 | 1,449 |
| DARTS-256 | 2,208 | 2,880 | 2,489 |
| DARTS-512 | 4,672 | 6,016 | 5,852 |

This version also improves hyperball sampling precision and fixes consistency issues between the reference and optimized implementations.

## Contents

- `Implementations/`: reference, optimized, and additional implementations.
- `Test_Vectors/`: known answer test vectors.
- `Algorithm specification (Version 2.0).pdf`: Version 2.0 algorithm specification.

See `Implementations/README.md` for implementation details.
