<!--
Copyright (c) 2026, NAKATA Maho
SPDX-License-Identifier: BSD-2-Clause
-->

# vMPLAPACK Routines

This document summarizes the public routine behavior for users of the prototype. The authoritative
specification remains `docs/SPEC.md`; if this document and the specification disagree, the
specification wins.

## Common Types

Verified routines return:

```cpp
template <class REAL>
struct Rmidrad {
    REAL mid;
    REAL rad;
    Rstatus status;
};
```

For `status == Rstatus::ok`, the real dot product is guaranteed to lie in:

```text
[mid - rad, mid + rad]
```

The guarantee is interval inclusion. It does not mean that `mid` is exact, and it does not require an
MPFR oracle at runtime.

`Rstatus` values:

```text
ok             finite enclosure constructed successfully
unbounded      finite inputs, but only the trivial infinite-radius enclosure is available
non_finite     at least one input is NaN or Inf
invalid_input  API boundary violation
```

The dot routines in this document use the same boundary rules:

```text
n < 0                                invalid_input
n == 0                               {0, 0, ok}
n > 0 and null x/y pointer            invalid_input
n > 0 and incx < 1 or incy < 1        invalid_input
finite inputs containing NaN or Inf   non_finite
```

Accurate routines such as `Rdot` have a lighter BLAS-style precondition contract; the explicit status
classification above is for verified routines.

## M9 Contract Validation Support

M9 does not add a new public numerical routine. It hardens the environment around the existing routines:

```text
configure-time guard   -ffast-math try_compile must fail through the umbrella __FAST_MATH__ #error
runtime sanity test    test_runtime_contract_m9
sanitizer build        VMPLAPACK_ENABLE_SANITIZERS=ON adds ASan+UBSan flags
CI matrix              native-only and MPFR-on, gcc and clang, Release and Debug
```

`test_runtime_contract_m9` is intentionally separate from the M0 environment smoke test. It checks the
ordinary accepted-flags build for runtime-directed rounding with constant folding blocked, x86 FTZ/DAZ
where MXCSR is available, and EFT invariants over normal, subnormal, cancellation, large/small-scale,
and signed-zero cases. MPFR-enabled builds run the same executable again at W=53 and W=512 to check the
MPFR rounding scope and representative EFT invariants.

These checks are not part of the library API, but they protect the assumptions used by `vRsum`,
`vRdot`, `vRdot_apriori`, and `vRresidual`.

## M10 Package Consumption Support

M10 does not add a numerical routine. It makes the existing routines consumable through an installed
CMake package:

```text
installed target        vmplapack::vmplapack
package config          vmplapackConfig.cmake
version policy          ExactVersion for 0.x releases
native-only package     no GMP/MPFR dependency lookup in find_package(vmplapack)
MPFR package            find_dependency(GMP) and find_dependency(MPFR), plus gmpfrxx_mkII headers
consumer smoke          tests/consumer builds and runs a verified-dot example
```

The installed `vmplapack::vmplapack` target carries the same strict FP usage requirements as the build
tree target. The consumer smoke project verifies this by inspecting `compile_commands.json`; it is not
assumed from successful compilation alone. The umbrella header remains the final defense against a
consumer overriding usage requirements with `-ffast-math`.

## M11 Examples And Benchmark Support

M11 does not add a new numerical kernel. It adds product-level demonstration and measurement drivers
around the existing public routines:

```text
example_m11_dot_driver         naive dot vs Rdot vs vRdot vs vRdot_apriori
example_m11_residual_worked    worked row-wise vRresidual example
benchmark_m11_core             CSV benchmark, smoke by default, full sweep opt-in
```

`example_m11_dot_driver` uses the M8 seeded high-condition generator and the MPFR oracle to display
condition estimates, oracle midpoint, midpoint errors, verified radii, statuses, and oracle inclusion.
The oracle is a reporting/checking tool only; the verified routines still compute their certificates
without oracle access.

`example_m11_residual_worked` constructs a two-row residual problem. One row has the cancellation shape
`[2^k, 1, -2^k] dot [1, 1, 1]`, so the example shows why row-wise verified residual boxes are more
meaningful than a displayed nearest-rounded residual alone.

`benchmark_m11_core` emits schema version `m11.1`. The schema is intentionally frozen for reuse by
M12+ tooling and includes timing, condition, status, enclosure, radius, compiler, platform, git, and
rounding-backend fields. Rows for routines without enclosures use explicit `null` for `enclosed`,
`radius`, and `relative_radius`; `unbounded` verified rows use `radius = +inf` and `enclosed = false`.
Smoke mode is intended for CI. Full mode adds extra random condition targets and orderings for manual
benchmark sweeps.

## M12a Verified Matvec/Matmul Reference

M12a adds the first verified linear-algebra layer status type:

```cpp
enum class VerificationStatus { Verified, Unverified, InvalidInput, Unsupported };
```

This status is separate from scalar `Rstatus`. Each output component is still an `Rmidrad<REAL>` with
its own scalar status. The top-level status summarizes the matrix operation:

```text
Verified       every requested component has Rstatus::ok
Unverified     arguments are well-formed, but at least one component is unbounded or non_finite
InvalidInput   dimension, pointer, stride, or leading-dimension contract violation
Unsupported    reserved for future LA modes not implemented by the current reference path
```

M12a point matvec, row-major storage:

```cpp
template <class REAL>
VerificationStatus vRgemv_point(std::ptrdiff_t m,
                                std::ptrdiff_t n,
                                const REAL* A,
                                std::ptrdiff_t lda,
                                const REAL* x,
                                std::ptrdiff_t incx,
                                Rmidrad<REAL>* out);
```

It computes enclosures for `out[i] ~= sum_j A[i*lda+j] * x[j*incx]`. Boundary rules:

```text
m < 0 or n < 0                         InvalidInput
m == 0                                  Verified, write nothing
m > 0 and out == nullptr                InvalidInput
n == 0                                  out[i] = {0, 0, ok}
n > 0 and A/x null, lda < n, incx < 1   InvalidInput
non-finite or overflowed component       component Rstatus, top-level Unverified
```

M12a point matmul, row-major storage:

```cpp
template <class REAL>
VerificationStatus vRgemm_point(std::ptrdiff_t m,
                                std::ptrdiff_t n,
                                std::ptrdiff_t k,
                                const REAL* A,
                                std::ptrdiff_t lda,
                                const REAL* B,
                                std::ptrdiff_t ldb,
                                Rmidrad<REAL>* C,
                                std::ptrdiff_t ldc);
```

It computes enclosures for `C[i*ldc+j] ~= sum_t A[i*lda+t] * B[t*ldb+j]`. Boundary rules mirror
`vRgemv_point`: negative dimensions are invalid, `m == 0 || n == 0` writes nothing and verifies,
`k == 0` writes exact zero boxes, and otherwise `C`, `A`, `B`, `ldc >= n`, `lda >= k`, and `ldb >= n`
are required.

M12a is the slow correctness baseline: each component can be certified by the existing directed
scalar `vRdot` path. It does not call BLAS/MPLAPACK and does not assume any external kernel honors
directed rounding.

## M12b Verified Matvec/Matmul Fast Path

The exported `vRgemv_point` and `vRgemm_point` signatures and boundary rules are unchanged. For each
component, the current implementation first tries the fast scalar a-priori certificate:

```text
mid = Rdot(component operands)                  nearest-rounding compensated Dot2 midpoint
rad = Dot2 a-priori radius using S_up           upward bound of sum |a_i*b_i|
```

This is the same contract as `vRdot_apriori`: the result is a rigorous `Rmidrad<REAL>` when the
component status is `ok`. The fast path uses nearest-rounding arithmetic for the midpoint and only
uses directed rounding for the scalar error-bound ingredients; it does not rely on BLAS/MPLAPACK
honoring a directed rounding mode.

If the a-priori path returns `unbounded`, the component falls back to the M12a directed `vRdot`
reference enclosure and keeps the better scalar status. This matters for cases such as huge exact
cancellation where `sum |a_i*b_i|` overflows in the a-priori bound but the directed enclosure of the
actual dot product is finite. `non_finite` and `invalid_input` are not hidden by the fallback.

`vRgemm_midrad` is still deferred; M12b covers point matrix inputs only.

## M13 Verified Solve Contract

M13 implements the public solve signatures below. Both overloads solve square point
systems and return a top-level `VerificationStatus`:

```cpp
template <class REAL>
VerificationStatus vRgesv(std::ptrdiff_t n,
                          const REAL* A,
                          std::ptrdiff_t lda,
                          const REAL* b,
                          std::ptrdiff_t incb,
                          Rmidrad<REAL>* x,
                          std::ptrdiff_t incx);

template <class REAL>
VerificationStatus vRgesv(std::ptrdiff_t n,
                          std::ptrdiff_t nrhs,
                          const REAL* A,
                          std::ptrdiff_t lda,
                          const REAL* B,
                          std::ptrdiff_t ldb,
                          Rmidrad<REAL>* X,
                          std::ptrdiff_t ldx);
```

Storage is row-major and all strides are measured in scalar elements. `A` is `n x n` with
`A[i*lda+j]` and `lda >= n`. In the vector overload, `b[i*incb]` is the input RHS and `x[i*incx]` is
the solution enclosure. In the matrix overload, `B[i*ldb+r]` is RHS column `r` and `X[i*ldx+r]` is
the matching solution enclosure, with `ldb >= nrhs` and `ldx >= nrhs`. Input and output ranges must
not overlap.

The certificate is normwise in the infinity norm. The routine may use any untrusted floating-point
inner engine to compute an approximate solution `x~`/`X~` and approximate inverse/preconditioner `R`.
It then verifies:

```text
alpha = upward_bound(||I - R*A||_inf)
beta  = upward_bound(||R*(b - A*x~)||_inf)   per RHS column
rad   = beta / (1 - alpha)
```

If `alpha < 1`, `A` is certified nonsingular and each output component for that RHS is
`{mid = x~_i, rad = rad, status = ok}`. For matrix RHS, `alpha` is shared across columns and `beta`
/ `rad` are per column. The residual `b - A*x~` must be computed accurately using the M12/M4 dot
infrastructure; a naive residual is not part of the M13 contract.

Current implementation notes:

- The approximate inverse/preconditioner `R` and approximate solution `x~`/`X~` are produced by a
  local nearest-rounding Gauss-Jordan engine with partial pivoting. This is an untrusted inner engine;
  vMPLAPACK does not link to MPLAPACK for the certificate.
- `alpha` is built from component enclosures of `R*A`: each component uses
  `vRdot_apriori_with_reference_fallback`, then contributes `|delta_ij - mid| + rad` under upward
  rounding to the row infinity-norm sum.
- Each residual component `b_i - A_i*x~` is computed from a verified dot enclosure for `A_i*x~` and
  an outward-rounded subtraction from `b_i`; the implementation does not use a naive residual.
- `beta` uses the interval residuals via the safe bound `|residual_j| <= |mid_j| + rad_j`, then sums
  `|R_ij|*(|mid_j| + rad_j)` under upward rounding. The resulting `rad = beta/(1-alpha)` is rounded
  upward and is shared by all components of that RHS column.

`example_m13_verified_solve` demonstrates all public M13 paths: vector RHS, matrix RHS,
near-singular-but-certified input, residual display for the returned midpoint, and a singular-looking
case that correctly reports `Unverified` with `unbounded` boxes. `example_m13_verified_solve_condition_sweep`
varies a 2x2 near-singular gap across each tier's precision boundary. `example_m13_verified_solve_random`
constructs a random strictly diagonally dominant integer system with a known exact solution.
`example_m13_verified_solve_hilbert` uses a Hilbert point matrix and a seeded integer solution vector
shape to display midpoint drift, residual boxes, and certificate radii for a harder classical family.
The random and Hilbert examples accept any `--n >= 1` so MPFR runs can push beyond the small native
smoke sizes; large orders are expensive because the dense solve and verification work scale cubically,
and the random example stores an `n` by `n` matrix. The random and Hilbert solve examples default the
MPFR tier to 512 bits and follow
`MPFRXX_DEFAULT_PRECISION_BITS` for lower-precision stress runs.

Boundary and status rules:

```text
n < 0, nrhs < 0, invalid pointers/strides/leading dimensions   InvalidInput
n == 0, or matrix overload nrhs == 0                            Verified, write nothing
non-finite A or RHS input                                       Unverified, non_finite boxes
alpha >= 1 or bound construction overflows/fails                Unverified, unbounded boxes
alpha < 1 and all requested boxes ok                            Verified
```

`Unverified` is only a failed certificate. It is never a claim that `A` is singular.

## `vRdot`

Signature:

```cpp
template <class REAL>
Rmidrad<REAL> vRdot(std::ptrdiff_t n,
                    const REAL* x,
                    std::ptrdiff_t incx,
                    const REAL* y,
                    std::ptrdiff_t incy);
```

`vRdot` constructs a rigorous directed-rounding enclosure for the dot product:

```text
z = exact sum_i x[i*incx] * y[i*incy]
```

It performs two plain product-and-add passes:

```text
sup = upward-rounded dot pass
inf = downward-rounded dot pass
```

Then it uses the accurate Dot2 routine as the midpoint:

```text
mid = Rdot(n, x, incx, y, incy)
return Rmake_midrad(inf, sup, mid)
```

If the accurate midpoint is not finite while the directed enclosure endpoints are finite,
`Rmake_midrad` falls back to the directed midpoint so the certificate remains valid.

Properties:

```text
soundness       direct interval enclosure from directed rounding
midpoint        accurate Dot2 result when finite
radius source   directed lower/upper pass width
oracle use      none at runtime
```

`vRdot` is the baseline verified dot. It is robust and simple, but its radius can be very wide for
strongly cancelling inputs. For example, `[2^k, 1, 2^k] dot [1, 1, -1]` has exact value `1`, but the
directed plain dot passes may produce a much wider enclosure around the cancellation.

## `vRdot_apriori`

Signature:

```cpp
template <class REAL>
Rmidrad<REAL> vRdot_apriori(std::ptrdiff_t n,
                            const REAL* x,
                            std::ptrdiff_t incx,
                            const REAL* y,
                            std::ptrdiff_t incy);
```

`vRdot_apriori` uses the same midpoint as `vRdot`:

```text
mid = Rdot(n, x, incx, y, incy)
```

The radius is not taken from directed lower/upper dot endpoints. It is a computable, outward-rounded
Ogita-Rump-Oishi Dot2 a-priori error bound:

```text
rad >= outward((u*abs(mid) + gamma_n^2*S_up + underflow_term) / (1 - u))
```

Definitions:

```text
p                 Rarith<REAL>::precision_bits()
u                 2^-p
gamma_n           n*u / (1 - n*u)
S_up              upward-rounded upper bound for sum_i abs(x_i*y_i)
underflow_term    5*n*eta, where eta = Rarith<REAL>::eta()
```

The exact ORO estimate is in terms of the unknown exact dot value `z`. Since `z` is not available at
runtime, the implementation uses `|z| <= |mid| + |mid - z|` and solves for the unknown error. This is
why the final computable bound has the division by `(1 - u)`.

`S_up` is computed by the named primitive:

```cpp
template <class REAL>
REAL Rsum_abs_dot_upper(std::ptrdiff_t n,
                        const REAL* x,
                        std::ptrdiff_t incx,
                        const REAL* y,
                        std::ptrdiff_t incy);
```

`Rsum_abs_dot_upper` runs in an upward-rounding scope and accumulates upward-rounded
`abs(x_i)*abs(y_i)`. It must not call accurate `Rdot`, because `Rdot` is an approximation and can
underestimate `sum_i abs(x_i*y_i)`.

Validity gate:

```text
n*u < 1
```

This is checked with integer arithmetic as:

```text
n < 2^p
```

If the gate fails, `vRdot_apriori` keeps the useful Dot2 midpoint but returns an infinite radius:

```text
{Rdot(...), +inf, unbounded}
```

MPFR underflow is checked explicitly. The MPFR tier does not assume that a wide exponent range makes
underflow impossible. If underflow is detected in the relevant bound path, the result is reported as
`unbounded`.

Properties:

```text
soundness       analytical Dot2 error bound, rounded outward
midpoint        accurate Dot2 result
radius source   ORO a-priori bound
oracle use      none at runtime
```

This routine is the M7 a-priori Dot2 estimate, not ORO Algorithm 5.8 `Dot2Err`. `Dot2Err` is a
different pure-floating-point error-bound routine with a different gate and is not implemented here.

## Choosing Between Them

Use `vRdot` when:

```text
you want the direct directed-rounding enclosure
you prefer the simpler verified mechanism
the problem size or platform behavior makes the a-priori gate/status less attractive
```

Use `vRdot_apriori` when:

```text
you want a Dot2-centered certificate
the input may have strong cancellation
you want a much smaller radius than the directed plain-dot enclosure can provide
```

In many high-condition cancellation examples, both routines are sound:

```text
vRdot          wider interval from directed plain dot endpoints
vRdot_apriori narrower interval from the Dot2 a-priori error bound
```

The narrower interval is not obtained by using an oracle. It follows from the Dot2 error theorem,
`S_up`, the length gate, the underflow policy, and outward rounding.

## Example Interpretation

For:

```text
[2^700, 1, 2^700] dot [1, 1, -1]
```

the exact value is:

```text
2^700 + 1 - 2^700 = 1
```

A naive nearest-rounded dot can return `0` because `2^700 + 1` rounds back to `2^700`. Dot2 can recover
the lost low-order contribution, so `Rdot` can return `1`. Therefore both `vRdot` and
`vRdot_apriori` may have:

```text
mid = 1
```

The important difference is the radius:

```text
vRdot.rad          may be very large because the directed plain-dot endpoints are far apart
vRdot_apriori.rad  can be tiny because it bounds the Dot2 error directly
```

When `vRdot_apriori` returns:

```text
mid = 1
rad = 5.26e-97
status = ok
```

the claim is:

```text
exact dot product is in [1 - 5.26e-97, 1 + 5.26e-97]
```

This statement is produced without using the MPFR oracle. The oracle appears only in tests and
diagnostic examples to check or explain the implementation.

## M8 Ill-Conditioned Generators

M8 adds generator routines in `vmplapack::gendot` for tests, diagnostics, and benchmark input
construction. These are not BLAS kernels and are not part of the core installed numerical API yet, but
they are product-level support routines for producing reproducible high-condition dot-product cases.
They require `VMPLAPACK_ENABLE_MPFR` because their metadata uses the MPFR oracle.

The basic generated dot container is:

```cpp
template <class REAL>
struct Rdot_case {
    std::vector<REAL> x;
    std::vector<REAL> y;
    REAL exact;
};
```

M8 metadata is carried by:

```cpp
enum class Rgenerator_status {
    ok,
    invalid_target,
    unachievable
};

enum class Rpermutation {
    generated,
    sorted,
    reversed,
    shuffled
};

template <class REAL>
struct Rgenerated_dot_case {
    Rdot_case<REAL> data;
    Rgenerator_status status;
    mpfrxx::mpfr_class target_cond_oro;
    mpfrxx::mpfr_class measured_cond_oro;
    mpfrxx::mpfr_class measured_cond_sum;
    std::uint64_t seed;
    int scale;
    Rpermutation permutation;
    const char* tier;
};
```

Status meaning:

```text
ok              case generated and measured condition metadata is populated
invalid_target  finite target condition is invalid under the ORO convention, such as target < 2
unachievable    target cannot be represented for the tier, length, exponent range, or construction
```

Condition-number convention:

```text
S        = sum_i abs(x_i*y_i)
cond_sum = S / abs(dot)
cond_oro = 2*S / abs(dot)
```

`cond_oro` is the ORO 2005 convention and is the target/benchmark axis. For nonzero dot products,
`cond_oro >= 2`. If the exact dot is zero and `S > 0`, both condition values are recorded as `+inf`.
This is a valid adversarial case, but not a finite log-scale target.

The measured condition fields are computed with the MPFR oracle interval and an upward high-precision
absolute-term pass. They are not computed with same-tier `Rdot` and not with a bare W-bit mpfr dot.

### Targeted Generator

```cpp
template <class REAL>
Rgenerated_dot_case<REAL> randomized_high_condition_power2(int target_cond_power2,
                                                           std::uint64_t seed,
                                                           Rpermutation permutation,
                                                           std::size_t pairs = 4U);
```

This constructs a multiset of cancelling `+/- 2^scale` terms plus one exact residual `1`. The target is

```text
target_cond_oro = 2^target_cond_power2
```

and `scale` is chosen so the measured ORO condition lands near the target. The generator uses the seed
to choose pair order in the generated multiset; `shuffled` additionally shuffles the full multiset.
`sorted` and `reversed` are ordering variants of the same values.

A successful finite target is accepted in log scale:

```text
abs(log10(measured_cond_oro) - log10(target_cond_oro)) <= 0.25
```

Tests and diagnostics can run this generator for several target condition powers per tier, such as low,
medium, and high finite targets, and compare each generated case with `vRdot` and `vRdot_apriori`. If
the target is below 2, the status is `invalid_target`. If the required scale is outside the tier's safe
exponent construction range, the status is `unachievable`.

### Deterministic Adversarial Families

```cpp
template <class REAL>
Rgenerated_dot_case<REAL> adversarial_exact_cancellation(int scale,
                                                         std::uint64_t seed,
                                                         Rpermutation permutation);

template <class REAL>
Rgenerated_dot_case<REAL> adversarial_heavy_cancellation(int scale,
                                                         std::uint64_t seed,
                                                         Rpermutation permutation);

template <class REAL>
Rgenerated_dot_case<REAL> adversarial_alternating_signs(int pairs,
                                                        int residual_exponent,
                                                        std::uint64_t seed,
                                                        Rpermutation permutation);

template <class REAL>
Rgenerated_dot_case<REAL> adversarial_huge_small_mixing(int large_exponent,
                                                        int small_exponent,
                                                        std::uint64_t seed,
                                                        Rpermutation permutation);
```

The families cover:

```text
exact cancellation        nontrivial S with exact dot 0, condition +inf
heavy cancellation        [2^scale, 1, -2^scale]-style cancellation
alternating signs         many +/-1 terms plus a small exact residual
huge/small scale mixing   large cancelling terms mixed with tiny cancelling terms and residual 1
```

All ordering variants must preserve the same multiset. They change evaluation order only; the exact dot
and condition convention do not change.

### Helper Routines

The main subordinate helpers are:

```cpp
template <class REAL>
REAL power_of_two(int exponent);

template <class REAL>
int max_power_exponent();

template <class REAL>
int min_power_exponent();

template <class REAL>
void measure_conditions(Rgenerated_dot_case<REAL>& result);

template <class REAL>
void apply_permutation(std::vector<Rterm<REAL>>& terms,
                       Rpermutation permutation,
                       std::uint64_t seed);
```

`power_of_two` is the exact-construction primitive: native tiers use in-range `std::ldexp`; MPFR uses
`mpfr_set_ui_2exp` at the current fixed precision `W`. `measure_conditions` is the oracle-based
metadata path. `apply_permutation` implements ordering variants without changing the multiset.

### How To Use The Metadata

Use `status` first. Only `status == ok` carries a generated case intended for numerical testing. Then
use:

```text
data.x, data.y          inputs to Rdot/vRdot/vRdot_apriori
measured_cond_oro       ORO condition number for plots and benchmark axes
measured_cond_sum       sum convention, useful for naive-dot first-order scaling
seed, scale             reproducibility metadata
permutation             generated/sorted/reversed/shuffled ordering label
```

The M8 tests check two theoretical scales:

```text
naive dot:   about gamma_n * cond_sum = 0.5 * gamma_n * cond_oro
Rdot/Dot2:   about u + 0.5 * gamma_n^2 * cond_oro
```

For every generated/adversarial case and every ordering, `vRdot` must enclose the MPFR oracle interval.
For every finite targeted high-condition random case across the selected target condition numbers,
`vRdot_apriori` is also compared with the oracle interval and must enclose it with `status == ok`. `Rdot`
is allowed to lose accuracy on high-condition inputs; verified inclusion is not allowed to fail.

## References

- T. Ogita, S. M. Rump, S. Oishi, "Accurate Sum and Dot Product", SIAM J. Sci. Comput. 26(6),
  1955-1988, 2005. DOI: 10.1137/030601818.
- `docs/SPEC.md`, M7: ORO a-priori accuracy bound.
