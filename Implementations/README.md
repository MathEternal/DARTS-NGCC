# DARTS Implementation

This project contains three different versions of the DARTS (lattice-based digital signature scheme) implementation: Reference Implementation, Optimized Implementation, and Additional Implementation. This project includes code derived from the HAETAE project, which is licensed under the MIT License. See Third Party License for details.

## Security status

The ordinary Reference, AVX2, and NEON trees use the SM3-based
`pseudohash`/`pseudoXOF` compatibility layer.  Its own API states that it
is provided only for implementation-correctness checks and does not guarantee
security.  Those trees must not be treated as production cryptographic
implementations.

The message-representative length is now independent of `CRHBYTES`:
DARTS-128/256 use 64-byte `mu`, while DARTS-512 uses 128-byte `mu`.  This is a
correctness-profile choice under the supplied auxiliary-hash rules, not a
claim that those auxiliary functions provide the corresponding collision
security.  The FIPS202 trees use SHAKE, but merely squeezing a longer SHAKE256
output does not raise SHAKE256 above its standardized 256-bit maximum
collision strength.  Key-generation seed-entropy and other design questions
also remain open.  Review this security status before using or evaluating any
parameter set.

## Project Structure

```
Implementations/
├── Reference_Implementation/         # Reference implementation
│   ├── DARTS128/
│   ├── DARTS256/
│   └── DARTS512/
├── Optimized_Implementation/         # Optimized implementation
│   ├── DARTS_AVX2/                   # AVX2 optimized
│   │   ├── DARTS128/
│   │   ├── DARTS256/
│   │   └── DARTS512/
│   └── DARTS_neon/                   # NEON optimized
│       ├── DARTS128/
│       ├── DARTS256/
│       └── DARTS512/
├── Additional_Implementation/        # Additional implementation (FIPS 202 standard)
│    └── DARTS_FIPS202/
│       ├── Implementations/
│       │   ├── Optimized_Implementation/
│       │   │   ├── DARTS128/
│       │   │   ├── DARTS256/
│       │   │   └── DARTS512/
│       │   └── Reference_Implementation/
│       │       ├── DARTS128/
│       │       ├── DARTS256/
│       │       └── DARTS512/
│       ├── KAT/                      # Known Answer Tests
│       └── helper_scripts/           # Helper scripts and tools
│
└── Third Party License/
     └── license.txt
```
## File Description

### Core Source Files

| File | Description |
|------|-------------|
| `sign.c/h` | Main logic for signature generation and verification |
| `poly.c/h` | Polynomial operations (addition, subtraction, multiplication, etc.) |
| `polyvec.c/h` | Polynomial vector operations |
| `polymat.c/h` | Polynomial matrix operations |
| `ntt.c/h` | Number Theoretic Transform (NTT) implementation |
| `fft.c/h` | Fast Fourier Transform (FFT) implementation |
| `sampler.c/h` | Sampler - sampling from discrete Gaussian distribution |
| `encoding.c/h` | rANS encoding/decoding |
| `packing.c/h` | Data packing and unpacking |
| `symmetric.c/h` / `symmetric-shake.c` | Symmetric cryptographic primitives (SHA3/SHAKE) |
| `fixpoint.c/h` | Fixed-point arithmetic |
| `polyfix.c/h` | Fixed-point polynomial operations |
| `fips202.c/h` | FIPS 202 standard implementation (Additional Implementation only) |
| `fips202x4.c/h` | FIPS 202 AVX2 optimized version (Additional Implementation only) |
| `consts.c/h` | Constant definitions and precomputed values (Optimized and Additional Implementation only) |
| `randombytes.c/h` | Random number generation |
| `reduce.c/h` / `reduce_avx.h` | Modular reduction operations |
| `drng.c/h` | Deterministic Random Number Generator |
| `auxfunc.c/h` | Auxiliary functions |
| `align.h` | Memory alignment definitions (Optimized Implementation) |
| `f1600x4.S` | SHA3-1600 assembly optimization (Optimized and Additional Implementation) |

### Header Files

| File | Description |
|------|-------------|
| `config.h` | Compilation configuration options |
| `params.h` | Security parameter definitions (different parameters for DARTS128/256/512) |

### Testing and Verification

| File | Description |
|------|-------------|
| `KAT_SIG.c` / `PQCgenKAT_sign.c` | Known Answer Tests - KAT generation program |
| `test/test_DARTS.c` | DARTS signature correctness test |
| `test/test_mul.c` | Polynomial multiplication test |
| `test/test_speed.c` | Performance benchmark test |
| `test/test_Sign.c` | Signature verification test (Additional Implementation) |
| `test/cpucycles.c/h` | Platform counter (ARM generic-timer ticks, x86 PMC/TSC counts, or fallback microseconds) |
| `test/speed_print.c/h` | Performance output formatting tool |

### Data Files

| Directory | Description |
|-----------|-------------|
| `precomputations_rANS/` | rANS precomputed data |
| `precomputations_rANS/prec_encoding.c` | rANS encoding precomputation code |
| `precomputations_rANS/rANS_*.txt` | rANS probability distribution tables |

### Other Files

| File | Description |
|------|-------------|
| `rans_byte.h` | Byte-level rANS implementation |
| `SIG_AlgorithmInstance.c/h` | Algorithm instance wrapper |
| `Makefile` | Build file |
| `README.txt` | Version-specific notes |
| `Third Party License/license.txt` | MIT License (Team HAETAE, 2023-2024) |

## Compilation and Running

### Reference Implementation (Reference_Implementation)

The reference implementation provides a basic, easy-to-understand implementation.

```bash
cd Reference_Implementation/DARTS128
make              # Compile all executables to bin/ directory
make test         # Generate test programs
make KAT          # Generate KAT program
make clean        # Clean build artifacts
```

Compiled executables are located in the `bin/` directory:
- `bin/KAT_SIG` - Known Answer Tests program
- `bin/test_DARTS` - DARTS signature test
- `bin/test_mul` - Polynomial multiplication test
- `bin/test_speed` - Performance benchmark test

### Recomputing the rANS model

The corrected statistics driver is shared by the three Reference parameter
sets. It uses valid keys, derives each signing seed as the production signer
does, applies the symmetric `(1-2b)c` Step-9 check, and retries candidates
that `pack_sig` cannot emit.

```bash
# Fast smoke test (repeat from DARTS128, DARTS256, and DARTS512)
make -B rans-stat CPPFLAGS=-DNUM_TESTS=1000
./bin/test_rANS

# Model-training run (example shown from Reference_Implementation/DARTS128)
make -B rans-stat \
  CPPFLAGS='-DNUM_TESTS=200000 -DKEY_ROTATION_INTERVAL=256 -DRANS_DRNG_SEED_BYTE=0x55'
./bin/test_rANS \
  > /tmp/darts_rans_train_128.txt \
  2> /tmp/darts_rans_train_128.progress

# Independent holdout, built only after installing the candidate model
make -B rans-stat \
  CPPFLAGS='-DNUM_TESTS=200000 -DKEY_ROTATION_INTERVAL=64 -DRANS_DRNG_SEED_BYTE=0xA5'
./bin/test_rANS \
  > /tmp/darts_rans_holdout_128.txt \
  2> /tmp/darts_rans_holdout_128.progress
```

Repeat both commands in `DARTS256` and `DARTS512`, changing the report names to
match the mode. Progress is written to standard error, leaving standard output
as the reproducible report. The `0x55` run, with one fresh valid key per 256
emitted signatures, is the only run used to derive `F_H` and `F_HB_Z1`. The
`0xA5` run deliberately uses a different deterministic stream and rotates keys
every 64 signatures; it is independent of the training stream and never
supplies the checked-in `F_*` arrays. If its result causes a cap or base change,
the pre-change report is tuning evidence rather than the final validation and
must not be presented as a report for the new configuration.

#### What the report counts

The report contains two related but different populations:

- The human-readable `z1` high-bit and `h` coefficient distributions near the
  start of the report contain exactly `NUM_TESTS` signatures that `pack_sig`
  emitted and the production verifier accepted. These emitted-only counts are
  useful for checking the final signing path, but they are not the input to
  frequency normalization.
- `COUNTS_H`, `COUNTS_HB_Z1`, the `F_*` lines, and all `PREPACK_*` size
  histograms cover every algorithm-accepted candidate whose complete
  coefficient vectors fit the production rANS alphabets, measured before
  `BASE_ENC_*` and `ENCODE_MAXSIZE` reject it. Consequently, a candidate that
  is representable but too long still trains the replacement model and appears
  in the uncensored size tail. The old table and size cap therefore cannot
  censor their replacement.

An algorithm-accepted candidate containing any coefficient outside the finite
production alphabet is reported through `UNMAPPED_*`, excluded as a whole from
`COUNTS_*`/`F_*`, and treated as a normal signing retry. Each representable
pre-pack candidate is directly encoded and decoded for the report. Each
emitted signature is separately verified through the production verifier.
`F_H` and `F_HB_Z1` must contain only positive integers and must each sum to
1024.

#### Fixed-point and holdout procedure

1. Run the three 200000-signature `0x55` training jobs and retain their complete
   reports and progress logs.
2. Validate and install the candidate frequencies and generated tables with
   `update_rans_tables.py --write`, initially without `--sync-reports`.
3. Rebuild and repeat the same `0x55` jobs with the installed tables. If either
   normalized frequency array changes, install it and repeat. The model is a
   fixed point only when a same-seed rerun reproduces both installed arrays
   exactly. A changed size cap or base also requires returning to this step,
   because it can change which candidates are needed to obtain `NUM_TESTS`
   emitted signatures.
4. With the fixed-point candidate installed, run the independent 200000-
   signature `0xA5` holdouts. Use `analyze_rans_reports.py` to validate the
   uncensored histograms, the current-cap tail, base windows, Step-9 rates, and
   final `b` balance. A cap selected from the empirical Q99.9 policy targets an
   observed tail no greater than 0.1%; it is an engineering criterion, not a
   proof about the unobserved distribution.
5. If holdout analysis changes `BASE_ENC_H`, `BASE_ENC_HB_Z1`,
   `ENCODE_MAXSIZE`, or the fixed signature size, synchronize those definitions
   first and repeat the fixed-point run and a fresh holdout under the final
   configuration. Do not reuse a report produced under a different cap as if
   it described the new configuration.

The analyzer is read-only and reads the current bases, cap, and fixed signature
portion from the Reference `params.h`. It recomputes histogram totals and
quantiles, verifies that the number of size rejections equals the observed tail
above the current cap, checks that component sizes fit their one-byte base
windows, and prints a candidate cap and `CRYPTO_SIGNATUREBYTES`:

```bash
python3 Additional_Implementation/DARTS_FIPS202/helper_scripts/analyze_rans_reports.py \
  --report 128=/tmp/darts_rans_holdout_128.txt \
  --report 256=/tmp/darts_rans_holdout_256.txt \
  --report 512=/tmp/darts_rans_holdout_512.txt
```

`analyze_rans_reports.py` does not alter tables or parameters. Its suggested
values require review and a new run under the selected configuration.

`update_rans_tables.py` has a different role: it strictly validates the three
model-input reports, recomputes the normalized frequencies from `COUNTS_*`,
regenerates all six encoder/decoder arrays per mode with the real C
precomputation program, and checks or updates all 15 `prec_encoding.c` and all
15 `encoding.c` copies:

```bash
python3 Additional_Implementation/DARTS_FIPS202/helper_scripts/update_rans_tables.py \
  --check \
  --report 128=/tmp/darts_rans_train_128.txt \
  --report 256=/tmp/darts_rans_train_256.txt \
  --report 512=/tmp/darts_rans_train_512.txt
```

The check is read-only and returns a nonzero status if any table is stale. By
default it rejects reports with fewer than 200000 successful samples, any
encode/decode failure, a wrong parameter tuple, a non-positive frequency, or a
frequency sum other than 1024. It does not choose bases or a size cap and does
not edit `params.h`, FIPS202 `api.h`, or KAT files.

After the final same-seed reports are fixed points, rerun the command with
`--write --sync-reports` to install the tables and copy each complete report to
all 15 `precomputations_rANS` directories (45 report replicas). Omitting
`--sync-reports` leaves those reports untouched. Adding it to `--check` compares
the 45 replicas without writing and is the required final consistency check.

#### Final synchronization checklist

A table, base, or cap change affects the wire format. Before declaring the
result final, synchronize and verify all of the following:

- rANS bases and `ENCODE_MAXSIZE` in every `params.h` copy;
- derived `CRYPTO_SIGNATUREBYTES` in `params.h` and the hard-coded value in all
  six FIPS202 `api.h` files;
- all 15 frequency sources, all 15 generated `encoding.c` files, and all 45
  report replicas through `update_rans_tables.py`;
- ordinary Reference, AVX2, and FIPS202 Reference/Optimized builds, direct
  encode/decode checks, signing self-tests, and cross-implementation checks;
- regenerated ordinary and FIPS202 KAT baselines, compared byte for byte across
  implementations. NEON execution must be performed on AArch64 hardware or a
  suitable cross-compiler/emulator environment.

The final numeric record below was produced by the fixed-point and independent
holdout runs described above.  The holdout column reports the uncensored tail
strictly above the installed cap, together with its empirical Q99.9 and
maximum.  Array normalization left both arrays unchanged for modes 128 and
512.  For mode 256, only zero-based `F_H` indices 4, 7, and 10 changed, from
`37,263,37` to `38,261,38`; `F_HB_Z1` was unchanged in all three modes.

| Mode | Fixed-point report SHA-256 | Holdout result | `BASE_ENC_H` | `BASE_ENC_HB_Z1` | `ENCODE_MAXSIZE` | `CRYPTO_SIGNATUREBYTES` | KAT regenerated |
|------|----------------------------|----------------|--------------|---------------------|------------------|-------------------------|-----------------|
| 128  | `f142579a7c9ee41d79057e0a29edea3a171654f95ec7c4e3af43297ae13b6bf5` | 118/200118 = 0.058965%; Q99.9=471; max=475 | 80 | 220 | 471 | 1449 | Yes |
| 256  | `0449dfaca4929af7ee9de007e4a95ab1dd73a09c6ed583c40942bd714471e359` | 80/200080 = 0.039984%; Q99.9=863; max=868 | 300 | 400 | 863 | 2489 | Yes |
| 512  | `5a7b6a5dc8b6f5bb3150dfad2fb01f472968258502de3d25f35c7d77d488032c` | 61/200061 = 0.030491%; Q99.9=2650; max=2653 | 1000 | 1500 | 2650 | 5852 | Yes |

### Optimized Implementation (Optimized_Implementation)

The optimized implementation provides AVX2 and NEON optimized versions to improve performance.

```bash
# AVX2 version
cd Optimized_Implementation/DARTS_AVX2/DARTS128
make              # Compile all executables to bin/ directory
make test         # Generate test programs
make KAT          # Generate KAT program
make clean        # Clean build artifacts
```

```bash
# NEON version
cd Optimized_Implementation/DARTS_neon/DARTS128
make              # Compile all executables to bin/ directory
make test         # Generate test programs
make KAT          # Generate KAT program
make clean        # Clean build artifacts
```

Compiled executables are located in the `bin/` directory, same as the reference implementation.

### Additional Implementation (Additional_Implementation - FIPS 202)

The additional implementation is based on the FIPS 202 standard and includes both optimized and reference implementations.

```bash
# Optimized version (FIPS 202)
cd Additional_Implementation/DARTS_FIPS202/Implementations/Optimized_Implementation/DARTS128
make              # Compile all executables to bin/ directory
make test         # Generate all test programs
make speed        # Generate performance test only
make clean        # Clean build artifacts
```

Compiled executables are located in the `bin/` directory:
- `bin/PQCgenKAT_sign` - KAT generation program
- `bin/test_speed` - Performance benchmark test
- `bin/test_Sign` - Signature verification test
- `bin/test_mul` - Polynomial multiplication test

## Parameter

DARTS provides three security levels of implementation:

| Version | Security Level | Parameter File |
|---------|----------------|----------------|
| DARTS128 | 128 bits | `params.h` |
| DARTS256 | 256 bits | `params.h` |
| DARTS512 | 512 bits | `params.h` |

The parameter differences for each version are defined in `params.h`, including polynomial modulus, coefficient ranges, etc.

## Examples

### Running KAT Tests

```bash
cd Reference_Implementation/DARTS128
make -B KAT
./bin/KAT_SIG
```

The ordinary Reference and AVX2 generators produce 10 deterministic vectors
per parameter set and verify every generated signature before writing it. The
checked-in baselines are in `../Test_Vectors`: files without a suffix belong
to the Reference implementation, while `_AVX2` files belong to AVX2. Compare
the generated file byte for byte with the matching baseline; a changed wire
format or rejection path requires a deliberately regenerated baseline, not a
relaxed comparison.

The FIPS202 `PQCgenKAT_sign` programs produce 100 vectors and call
`crypto_sign_open` for every vector. Their common Reference/Optimized
baselines are in `Additional_Implementation/DARTS_FIPS202/KAT`. In addition
to each implementation's self-test, Reference and optimized signatures must
be verified in both directions. Their deterministic public and secret keys
must match for the same request file.

Run each generator twice and use `cmp` to check reproducibility before
updating a baseline. NEON KAT execution requires an ARMv8/AArch64 machine (or
an AArch64 cross-compiler and emulator); a successful x86 build is not a
substitute for executing the NEON binary.

### Running Performance Benchmarks

```bash
cd Reference_Implementation/DARTS128
make
./bin/test_speed
```

### Verifying Signature Correctness

```bash
cd Reference_Implementation/DARTS128
make test
./bin/test_DARTS
```
