<!--
Copyright (c) 2026, NAKATA Maho
SPDX-License-Identifier: BSD-2-Clause
-->

# vMPLAPACK

vMPLAPACK is a standalone C++17 prototype of accurate and verified `sum` and `dot` kernels in the MPLAPACK style. It is not part of MPLAPACK and does not depend on the MPLAPACK source tree.

The public routines are templated on `REAL` and share one generic implementation across three co-equal scalar tiers:

```text
REAL = float
REAL = double
REAL = mpfrxx::mpfr_class
```

Tier-specific behavior is isolated in `vmplapack::Rarith<REAL>`. The MPFR tier uses gmpfrxx_mkII at a fixed working precision `W`; it is a first-class scalar tier, and high-precision MPFR is also used by the test oracle.

## Accurate vs Verified

Accurate routines return a bare `REAL` approximation:

```cpp
Rsum(n, x, incx)
Rdot(n, x, incx, y, incy)
```

They use compensated Sum2/Dot2-style algorithms. They are better approximations, but they do not carry a certificate.

Verified routines return a midpoint-radius enclosure plus a status:

```cpp
vRsum(n, x, incx)
vRdot(n, x, incx, y, incy)
vRresidual(m, n, A, lda, x, b, out)
vRgemv_point(m, n, A, lda, x, incx, out)
vRgemm_point(m, n, k, A, lda, B, ldb, C, ldc)
vRgesv(n, A, lda, b, incb, x, incx)
vRgesv(n, nrhs, A, lda, B, ldb, X, ldx)
```

For scalar `Rmidrad` components with `status == ok`, the true finite result is guaranteed to be in:

```text
[mid - rad, mid + rad]
```

For `status == unbounded`, the certificate is still valid but useless: `rad == +inf`, representing `[-inf,+inf]`. This is used when finite inputs overflow during bound construction. For `status == non_finite`, an input was NaN or Inf, so no finite-real certificate is claimed. For `status == invalid_input`, the arguments violate the verified API boundary rules.

Linear-algebra routines return `VerificationStatus` at the operation level: `Verified` means every
component box is `ok`; `Unverified` means at least one component is `unbounded` or `non_finite`;
`InvalidInput` means a boundary-rule violation; `Unsupported` is reserved for future LA modes.
`vRgemv_point` and `vRgemm_point` use the M12b a-priori component enclosure first and fall back to the
M12a directed component enclosure when the a-priori bound is unbounded. `vRgesv` returns a normwise
verified enclosure for square point systems when it can certify `||I - R*A||_inf < 1`.

Verification means interval inclusion, not closeness to an MPFR point. Tests compare verified intervals against an MPFR oracle interval `[ref_lo, ref_hi]`.

See `ROUTINES.md` for the user-facing behavior of `vRdot` and `vRdot_apriori`, including how to
interpret `mid`, `rad`, and `status`.

## Boundary Rules

Accurate routines follow the bare BLAS-style contract:

```text
n <= 0       -> return zero and inspect no pointer
n > 0        -> pointers valid and strides >= 1 are preconditions
```

Verified routines classify invalid inputs explicitly:

```text
vRsum/vRdot:
  n < 0       -> invalid_input
  n == 0      -> {0, 0, ok}
  null pointer or non-positive stride with n > 0 -> invalid_input

vRresidual:
  m < 0 or n < 0 -> invalid_input
  m == 0         -> ok, write nothing
  m > 0, n == 0  -> out[i] = {b[i], 0, ok} for finite b[i]
  m > 0, n > 0   -> A, x, b, out valid and lda >= n required

vRgesv:
  negative dimensions, bad pointers, bad strides, or bad leading dimensions -> InvalidInput
  n == 0, or matrix-RHS nrhs == 0 -> Verified, write nothing
  non-finite A or RHS input -> Unverified with non_finite boxes
  failed nonsingularity/bound certificate -> Unverified with unbounded boxes
```

`vRresidual` computes `r = b - A*x` row by row and returns the worst `Rstatus` by `Rstatus_rank`.

## Floating-Point Contract

The verified algorithms depend on runtime directed rounding. Native builds must preserve IEEE-754 semantics and must not contract or reassociate floating-point expressions. Consumers should link the CMake target rather than only adding include paths:

```cmake
target_link_libraries(your_target PRIVATE vmplapack::vmplapack)
```

The target propagates the required flags on supported compilers, including:

```text
-fno-fast-math
-frounding-math
-ffp-contract=off
```

The umbrella header rejects fast-math builds with `__FAST_MATH__`, and it requires `FLT_EVAL_METHOD == 0`. Runtime tests also check that rounding scopes restore the prior mode, subnormals are alive, and directed rounding affects ordinary addition and multiplication.

Accurate/EFT kernels require ambient round-to-nearest. Do not call them from inside a directed-rounding scope. Verified passes use plain product then add, not FMA contraction.

## MPFR Fixed Precision

When `VMPLAPACK_ENABLE_MPFR=ON`, `mpfrxx::mpfr_class` is enabled as a first-class tier. The working precision `W` is fixed for a test process or run. Set it before constructing any MPFR inputs:

```cpp
mpfrxx::set_default_precision_bits(512);
```

Existing `mpfr_class` objects keep their original precision, so do not mix objects created before and after changing `W`. Tests run MPFR at `W=53` and `W=512` as separate CTest cases.

## Build and Test

Full build with MPFR:

```bash
cmake -S . -B build -DVMPLAPACK_ENABLE_MPFR=ON -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Native-only build, with no GMP/MPFR requirement:

```bash
cmake -S . -B build-native -DVMPLAPACK_ENABLE_MPFR=OFF -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build-native -j
ctest --test-dir build-native --output-on-failure
```

M9 contract and sanitizer checks:

```bash
# The default configure runs a try_compile that must reject -ffast-math via the umbrella header.
cmake -S . -B build -DVMPLAPACK_ENABLE_FAST_MATH_GUARD_TEST=ON

# Debug ASan+UBSan build.
cmake -S . -B build-sanitize -DVMPLAPACK_ENABLE_SANITIZERS=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build build-sanitize -j
ASAN_OPTIONS=detect_leaks=0 ctest --test-dir build-sanitize --output-on-failure
```

Optional Dekker TwoProduct path:

```bash
cmake -S . -B build-dekker -DVMPLAPACK_USE_DEKKER_TWOPRODUCT=ON
```

If gmpfrxx_mkII is not in the default probe location, pass:

```bash
-DGMPFRXX_MKII_INCLUDE_DIR=/path/to/gmpfrxx_mkII/include
```

## Testing vRgemm_point

`vRgemm_point` is tested in the M12 suites. These focused tests are registered in MPFR-enabled
builds because they compare each component enclosure against an MPFR oracle interval; the same test
executables also exercise the native `float` and `double` tiers.

Build and run only the M12 `vRgemv_point` / `vRgemm_point` tests:

```bash
cmake -S . -B build -DVMPLAPACK_ENABLE_MPFR=ON -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j
ctest --test-dir build -R 'test_M12a_reference|test_M12b_fast_matmul' --output-on-failure
```

`test_M12a_reference` checks the reference directed-rounding enclosure, row-major leading-dimension
layout, boundary rules, `k == 0`, invalid inputs, and oracle inclusion. `test_M12b_fast_matmul`
checks the fast a-priori component enclosure, cancellation-heavy components, fallback behavior, and
that the fast path is tighter than the directed reference where expected.

The M12/M13 linear-algebra examples are useful for inspecting the actual boxes:

```bash
cmake --build build -j --target \
  example_m12_verified_gemm \
  example_m12_verified_gemm_midpoint_error \
  example_m12_verified_gemm_highcondition \
  example_m12_verified_gemm_medium_diff \
  example_m13_verified_solve \
  example_m13_verified_solve_condition_sweep \
  example_m13_verified_solve_random \
  example_m13_verified_solve_hilbert \
  example_m14_verified_inverse \
  example_m14_verified_inverse_highcondition

MPFRXX_DEFAULT_PRECISION_BITS=512 ./build/example_m12_verified_gemm
MPFRXX_DEFAULT_PRECISION_BITS=512 ./build/example_m12_verified_gemm_midpoint_error --m 4 --n 4 --k 10 --seed 17
MPFRXX_DEFAULT_PRECISION_BITS=512 ./build/example_m12_verified_gemm_highcondition --m 5 --n 5 --k 13 --seed 123
MPFRXX_DEFAULT_PRECISION_BITS=512 ./build/example_m12_verified_gemm_medium_diff --m 5 --n 5 --k 13 --seed 123
MPFRXX_DEFAULT_PRECISION_BITS=512 ./build/example_m13_verified_solve
MPFRXX_DEFAULT_PRECISION_BITS=512 ./build/example_m13_verified_solve_condition_sweep
MPFRXX_DEFAULT_PRECISION_BITS=512 ./build/example_m13_verified_solve_random --n 6 --seed 123
MPFRXX_DEFAULT_PRECISION_BITS=512 ./build/example_m13_verified_solve_hilbert --n 30
MPFRXX_DEFAULT_PRECISION_BITS=512 ./build/example_m14_verified_inverse
MPFRXX_DEFAULT_PRECISION_BITS=512 ./build/example_m14_verified_inverse_highcondition --n 20 --float-cond 12 --double-cond 40 --mpfr-cond 180
```

For quiet smoke runs, pass `--no-matrices` to the random/medium-difference examples:

```bash
MPFRXX_DEFAULT_PRECISION_BITS=512 ./build/example_m12_verified_gemm_midpoint_error --no-matrices
MPFRXX_DEFAULT_PRECISION_BITS=512 ./build/example_m12_verified_gemm_highcondition --no-matrices
MPFRXX_DEFAULT_PRECISION_BITS=512 ./build/example_m12_verified_gemm_medium_diff --no-matrices
MPFRXX_DEFAULT_PRECISION_BITS=512 ./build/example_m13_verified_solve_random --n 6 --seed 123 --no-matrices
MPFRXX_DEFAULT_PRECISION_BITS=128 ./build/example_m13_verified_solve_random --n 32 --seed 123 --no-matrices
MPFRXX_DEFAULT_PRECISION_BITS=512 ./build/example_m13_verified_solve_hilbert --n 30 --no-matrices
MPFRXX_DEFAULT_PRECISION_BITS=128 ./build/example_m13_verified_solve_hilbert --n 20 --no-matrices
MPFRXX_DEFAULT_PRECISION_BITS=512 ./build/example_m14_verified_inverse
MPFRXX_DEFAULT_PRECISION_BITS=512 ./build/example_m14_verified_inverse_highcondition --n 20 --float-cond 12 --double-cond 40 --mpfr-cond 180
```

`example_m12_verified_gemm_midpoint_error` is intentionally not a high-condition generator. It uses
`--m`, `--n`, `--k`, `--seed`, and `--no-matrices`, then fills both point matrices with independent
seeded signed rational entries. Each entry is `sign * n / d`, with `n` in `[1, 127]` and odd `d` in
`[65, 255]`, rounded in the active tier. This gives natural-looking small matrices and keeps
`C.diff = C.mid - oracle midpoint` at ordinary rounding-error scale.

`example_m12_verified_gemm_highcondition` is the condition-targeting example. It accepts `--cond`,
`--float-cond`, `--double-cond`, and `--mpfr-cond`. For each cancellation pair, `A` receives the
same high-magnitude dyadic random value in both pair positions, while `B` receives opposite signed
dyadic random values in the same two positions. The pair contribution is therefore exactly zero in
the mathematical dot product, but its absolute term sum is large. One signed rational residual term
per row/column keeps the exact product nonzero and less hand-picked, so the printed measured
condition is high while the oracle interval remains meaningful.

Read the output as follows: `return status = Verified` means every component has `Rstatus::ok`;
`all oracle intervals covered = true` is the MPFR-oracle inclusion check; `result C.rad` is the
verified certificate radius; and `result C.diff` is only the displayed center error
`C.mid - oracle C midpoint`, not the certificate. The random GEMM diagnostics also print
`max |C.mid - oracle midpoint| / u`, `max radius / u`, and worst-component `lower-mid` /
`upper-mid` offsets, so MPFR@512 intervals remain visible even when decimal endpoint printing makes
`[mid-rad, mid+rad]` look degenerate.

## Install and Consume

M10 installs an exported CMake package. Native-only installs do not require GMP/MPFR at consumer
configure time; MPFR-enabled installs call `find_dependency(GMP)` and `find_dependency(MPFR)` and need
gmpfrxx_mkII headers available through `GMPFRXX_MKII_INCLUDE_DIR` if the build-time path is not valid.

```bash
cmake --install build --prefix build/install
cmake -S tests/consumer -B build/consumer \
  -DCMAKE_PREFIX_PATH=$PWD/build/install \
  -DGMPFRXX_MKII_INCLUDE_DIR=/path/to/gmpfrxx_mkII/include \
  -DVMPLAPACK_CONSUMER_EXPECT_MPFR=ON
cmake --build build/consumer -j
ctest --test-dir build/consumer --output-on-failure
```

The consumer smoke project also checks `compile_commands.json` to verify that the strict FP flags are
inherited through `vmplapack::vmplapack`. A consumer that adds `-ffast-math` is rejected by the umbrella
header guard.

## Examples

Examples are built when `VMPLAPACK_ENABLE_EXAMPLES=ON`. Tier-complete examples exercise `float`, `double`, and `mpfrxx::mpfr_class` at `W=512`.

Useful examples include:

```text
example_vRdot
example_vRresidual
example_m4_accurate_dot
example_m4_accurate_sum
example_m7_apriori_bound
example_m11_dot_driver
example_m11_residual_worked
example_m12_verified_gemm
example_m12_verified_gemm_midpoint_error
example_m12_verified_gemm_highcondition
example_m12_verified_gemm_medium_diff
example_m13_verified_solve
example_m13_verified_solve_condition_sweep
example_m13_verified_solve_random
example_m13_verified_solve_hilbert
example_m14_verified_inverse
example_m14_verified_inverse_highcondition
example_m3_oracle_generators
```

`example_m11_dot_driver` compares naive dot, `Rdot`, `vRdot`, and `vRdot_apriori` on one seeded
high-condition dot per tier. It prints oracle midpoint, condition estimates, errors, radii, statuses,
and whether verified boxes enclose the oracle interval. `example_m11_residual_worked` shows a small
`vRresidual` problem with cancellation in one row and prints the row-wise residual boxes.
`example_m13_verified_solve` tours `vRgesv`: exact vector solve, matrix RHS solve, near-singular
certified solve, residual of the returned midpoint, and an `Unverified` certificate-failure case.
`example_m13_verified_solve_condition_sweep` sweeps the gap in a 2x2 near-singular system to show
where the stored point matrix stops being certifiable. `example_m13_verified_solve_random` builds a
random strictly diagonally dominant integer system with known exact solution and reports coverage,
midpoint error, and radius scale. `example_m13_verified_solve_hilbert` uses Hilbert matrices to show
larger midpoint drift and radii on a classical ill-conditioned family. Both examples accept any
`--n >= 1`; large `n` is intentionally allowed for MPFR experiments, but the dense solve and
certificate work scale cubically and the random example stores an `n` by `n` matrix. The random and
Hilbert solve examples default the MPFR tier to 512 bits and follow
`MPFRXX_DEFAULT_PRECISION_BITS`, so
`MPFRXX_DEFAULT_PRECISION_BITS=128 ./build/example_m13_verified_solve_random --n 32 --seed 123 --no-matrices`
and
`MPFRXX_DEFAULT_PRECISION_BITS=128 ./build/example_m13_verified_solve_hilbert --n 20 --no-matrices`
are useful low-precision stress runs. `example_m14_verified_inverse` prints verified inverse boxes for a
regular matrix and a near-singular 2x2 case, plus a verified bound for `A * inverse.mid - I`.
`example_m14_verified_inverse_highcondition` builds larger dense symmetric point matrices as
`Q * diag(1 ... 2^-p) * Q^T`, where `Q` is a deterministic product of signed Householder reflectors.
The command-line `--n`, `--cond`, `--float-cond`, `--double-cond`, and `--mpfr-cond` options control
the size and target spectral condition scale `2^p`. This is a generator target, not a certified
condition routine; M14b owns certified condition bounds.

## Benchmarks

M11 adds `benchmark_m11_core`, built when `VMPLAPACK_ENABLE_BENCHMARKS=ON` and MPFR support is enabled.
The default CTest path runs only the light smoke benchmark:

```bash
./build/benchmark_m11_core --smoke
./build/benchmark_m11_core --full
```

The benchmark emits CSV to stdout with frozen schema version `m11.1`:

```text
schema_version,routine,tier,precision_bits,n,generator,seed,target_cond,realized_cond_oro,
realized_cond_sum,status,enclosed,elapsed_ns_total,time_ns_per_item,work_items,repetitions,
statistic,mid_error_abs,mid_error_rel,radius,relative_radius,compiler_id,compiler_version,
build_type,fp_contract_flags,cpu,os,git_sha,rounding_backend
```

Rows without an enclosure, such as naive dot and `Rdot`, use `enclosed = null`, `radius = null`, and
`relative_radius = null`. Verified rows use `status`, `enclosed`, `radius`, and `relative_radius`;
`unbounded` rows use `radius = +inf` and `enclosed = false`. The cases include deterministic M3 data
and random/adversarial M8 data, so grouping by `tier`, `routine`, and `realized_cond_oro` gives
success-rate-vs-condition and radius-vs-condition views.

## References

- T. Ogita, S. M. Rump, S. Oishi, "Accurate Sum and Dot Product", SIAM J. Sci. Comput. 26(6), 1955-1988, 2005. DOI: 10.1137/030601818.
- S. M. Rump, "Verification methods: Rigorous results using floating-point arithmetic", Acta Numerica 19, 287-449, 2010. DOI: 10.1017/S096249291000005X.
- N. J. Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., SIAM, 2002. ISBN: 0-89871-521-0 / 978-0-89871-521-7. DOI: 10.1137/1.9780898718027.
