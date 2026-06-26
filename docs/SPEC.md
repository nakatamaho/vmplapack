<!--
Copyright (c) 2026, NAKATA Maho
SPDX-License-Identifier: BSD-2-Clause
-->

# Task: vMPLAPACK ― verified / accurate `sum` and `dot` kernels (standalone prototype)

You are an autonomous coding agent. Build **vMPLAPACK**: a small, self-contained C++17
library of **accurate** and **verified** reductions (`sum`, `dot`), in the MPLAPACK idiom,
templated on a real scalar type `REAL`. Three scalar tiers are **co-equal first-class
citizens**, unified by a single `Rarith<REAL>` policy:

```text
REAL = float                  // hardware IEEE-754 binary32          (u = 2^-24)
REAL = double                 // hardware IEEE-754 binary64          (u = 2^-53)
REAL = mpfrxx::mpfr_class      // MPFR via gmpfrxx_mkII, FIXED working precision W
```

All three tiers have a correctly-rounded FMA and runtime directed rounding, so the *same*
generic algorithms run on every tier. This is a sandbox prototype ― **not** an MPLAPACK patch.

This spec has been hardened against subtle soundness traps (operand aliasing, overflow/underflow
of bounds, NaN-vs-Inf classification, empty/invalid inputs, MPFR fixed-precision semantics, oracle
rounding, condition-aware testing). **These hardening rules are load-bearing for the word
"verified" ― do not relax them.**

Implement **milestone by milestone (M0 → M12)**, one commit/PR per milestone; do **not** start a
milestone until the previous one's acceptance tests pass. If you find a reason to deviate from the
signatures, `Rmidrad`/`Rstatus`, the `Rarith<REAL>` interface, the preconditions (§6), the
boundary rules (§8), or the rounding contract (§9), **stop and report** instead of silently
changing them.

---

## 0. Scope

In scope: `REAL`-generic EFTs (`Rtwosum`, `Rfasttwosum`, `Rtwoproduct`, `Rsplit`); accurate
compensated reductions (`Rsum`, `Rdot`); verified reductions returning a rigorous midrad enclosure
(`vRsum`, `vRdot`); a verified residual example (`vRresidual`); the `Rarith<REAL>` policy unifying
the three tiers; an MPFR naive oracle (tests only); a CMake build with a strict rounding/FP
contract; explicit handling of overflow / underflow / non-finite / empty inputs via status.

Out of scope: full MPBLAS/MPLAPACK; LU/QR/eigen/SVD or any verified solver; integration into the
real MPLAPACK tree; `long double` (not an expression scalar leaf in gmpfrxx_mkII; exclude it);
variable/per-object mpfr precision (the mpfr tier runs at a single FIXED `W`); performance tuning,
SIMD, threading, GPU.

---

## 1. Two concepts, kept distinct

- **Accurate** (`Rsum`, `Rdot`): a *better approximation* ( twice-working-precision via
  compensation). **No certificate.** Returns a bare `REAL`.
- **Verified** (`vRsum`, `vRdot`, `vRresidual`): an approximation **plus a rigorous, guaranteed
  bound** as a midpointradius (midrad) enclosure `true ∈ [mid - rad, mid + rad]`, `rad >= 0`,
  together with a `Rstatus`. "Close to the MPFR value" is **not** verification; interval
  **inclusion** is.

---

## 2. Architecture: one generic algorithm, three tiers via `Rarith<REAL>`

Do **not** write per-type kernels or a virtual per-scalar-op backend. Kernels are
`template <class REAL>`; **all** tier-specific behavior lives in the traits policy `Rarith<REAL>`
(the parallel of MPLAPACK's `get_arithmetic_params<REAL>()`).

`Rarith<REAL>` must provide:

```text
precision_bits()       -> long   // significand bits p; MUST be a SIGNED type (see §2.2 note)
unit_roundoff()        -> REAL   // u = 2^-p
zero()  one()  half()  -> REAL   // exact constants at the working precision
infinity()             -> REAL   // +inf at the working precision (for unbounded radius)
is_finite(x)           -> bool   // false for NaN AND Inf (the only finiteness query needed)
abs(x)                 -> REAL
fma(a, b, c)           -> REAL    // CORRECTLY-ROUNDED a*b + c
struct round_up   { ... };        // RAII directed rounding toward +inf; NON-COPYABLE, NON-MOVABLE
struct round_down { ... };        // RAII directed rounding toward -inf; NON-COPYABLE, NON-MOVABLE
```

`is_finite` is the **only** finiteness predicate; the design never needs to distinguish NaN from
Inf (see §6.3). `round_up`/`round_down` must delete copy and move (a copied scope would restore
twice). The accurate path runs in the ambient default mode, required to be round-to-nearest (§9.3).

### 2.1 Native tiers: `Rarith<float>` and `Rarith<double>`

Identical structure; differ only in constants and the overloaded `std::fma`/`std::isfinite`/
`std::fabs`. `fma` → `std::fma` (`fmaf` for `float`; correctly rounded; explicit call, immune to
`-ffp-contract`). `round_up`/`round_down` → RAII over `std::fesetround(FE_UPWARD/FE_DOWNWARD)`,
restoring the prior mode in the destructor; copy/move deleted. `float`: p=24, u=`0x1p-24f`;
`double`: p=53, u=`0x1p-53`. `half/one`: 0.5/1.0 in the tier type. `infinity`:
`std::numeric_limits<REAL>::infinity()`. `is_finite`→`std::isfinite`; `abs`→`std::fabs`.
`precision_bits()` returns a signed `long`.

### 2.2 `Rarith<mpfrxx::mpfr_class>` (verified against the gmpfrxx_mkII source; still probe ― §9.3 / M2)

Header `#include <gmpfrxx_mkII.h>` (or MPFR-only `#include <mpfrxx_mkII.h>`). Namespace `mpfrxx`;
type `mpfrxx::mpfr_class`. Working precision `W` is FIXED (§10). API to use:

- `fma` → free function `mpfrxx::fma(a, b, c)` (wraps `mpfr_fma`, correctly rounded).
- `round_up`/`round_down` → wrap `mpfrxx::rounding_mode_scope` (RAII; sets the thread's MPFR default
  rounding mode and restores it; already non-copyable):
  `struct round_up { mpfrxx::rounding_mode_scope s_{MPFR_RNDU}; };` and `MPFR_RNDD` for `round_down`.
- `precision_bits()` → `static_cast<long>(mpfrxx::default_prec())` (== `W`); **signed**.
- `unit_roundoff()`: build `2^-W` with a SIGNED exponent ― never write `-precision_bits()` if it
  could be unsigned:
  ```cpp
  mpfrxx::mpfr_class u = mpfrxx::mpfr_class::with_precision(precision_bits());
  mpfr_exp_t e = -static_cast<mpfr_exp_t>(precision_bits());   // signed before negation
  mpfr_set_ui_2exp(u.mpfr_data(), 1, e, MPFR_RNDN);
  ```
  Raw `mpfr_t` accessor: `mpfr_class::mpfr_data()` (also `get_mpfr_t()`).
- `zero()/one()/half()` → `mpfrxx::mpfr_class::with_precision(W, 0.0/1.0/0.5)` (all exact).
- `infinity()` → an `mpfr_class` at precision `W` with `mpfr_set_inf(_, +1)`.
- `is_finite(x)` → `mpfr_number_p(x.mpfr_data()) != 0`; `abs` → `mpfrxx::abs`.

Gate the whole specialization behind `#ifdef VMPLAPACK_ENABLE_MPFR` so a native-only build needs no
GMP/MPFR. When ON, the library links MPFR (mpfr is first-class here, not a test-only oracle).

### 2.3 Expression-template + materialization rule (ALL tiers)

gmpfrxx_mkII operators return **expression nodes**, not `mpfr_class`. **Materialize each elementary
operation into a named `REAL`** (one rounded op per statement). Never bind an arithmetic result with
`auto`; expand multi-op expressions into single-op assignments (see EFTs and Dot2). Free for native;
for mpfr it forces per-op rounding under the active mode.

### 2.4 Naming rules (MPLAPACK-faithful; no concrete type/library in names)

`R`-prefix = "real, templated over `REAL`" (family marker, not a concrete type). Accurate routines
keep MPLAPACK names (`Rsum`, `Rdot`); verified routines carry a leading `v` (`vRsum`, `vRdot`,
`vRresidual`). No `float`/`double`/`mpfr`/`gmp`/`native` in public names. Namespace `vmplapack`;
oracle (tests) `vmplapack::oracle`.

---

## 3. Repository layout (MPLAPACK-style: one routine per `.cpp`)

```text
vmplapack/
  CMakeLists.txt
  cmake/RoundingControl.cmake        # rounding/FP flag contract (§9), exposed PUBLIC/INTERFACE
  include/vmplapack/
    vmplapack.h                      # umbrella; contains the __FAST_MATH__ guard (§9.1)
    vmpblas.h                        # public decls + extern template decls for the first-class tiers
    vmplapack_utils.h                # Rstatus, Rstatus_rank, Rmidrad<REAL>, midrad helpers
    vmplapack_arith.h                # Rarith<REAL>
    vmplapack_eft.h                  # Rtwosum / Rfasttwosum / Rtwoproduct / Rsplit (template defs)
  src/
    Rsum.cpp  Rdot.cpp               # explicit instantiation: float, double, (mpfr if enabled)
    vRsum.cpp vRdot.cpp vRresidual.cpp
  tests/
    Rdot_oracle.h                    # naive high-precision mpfr oracle, returns an INTERVAL (§10)
    Rgendot.h                        # deterministic ill-conditioned generators (§11.4)
    test_Rdot.cpp  test_vRdot.cpp  test_vRresidual.cpp
    test_Rarith_environment.cpp      # rounding/FP contract: static + runtime checks
  examples/
    example_vRdot.cpp  example_vRresidual.cpp
  README.md
```

Template definitions live in the public headers (generic use with any conforming `REAL`). For the
first-class tiers, also provide **explicit instantiation** in `src/*.cpp` and **`extern template`**
declarations in `vmpblas.h` (guard the mpfr instantiation with `VMPLAPACK_ENABLE_MPFR`). No
dependency on the real MPLAPACK source tree.

---

## 4. Result types (`include/vmplapack/vmplapack_utils.h`)

```cpp
#pragma once
#include <cstddef>

namespace vmplapack {

enum class Rstatus {
    ok,             // finite mid and finite rad; inclusion certificate valid
    unbounded,      // inclusion certificate valid but rad is +inf (bound/endpoint/radius overflow)
    non_finite,     // a NaN/Inf INPUT; no certificate is claimed
    invalid_input   // bad length/stride/pointer
};

// Aggregation order for multi-row results: ok < unbounded < non_finite < invalid_input.
inline int Rstatus_rank(Rstatus s) {
    switch (s) {
        case Rstatus::ok:            return 0;
        case Rstatus::unbounded:     return 1;
        case Rstatus::non_finite:    return 2;
        case Rstatus::invalid_input: return 3;
    }
    return 3;
}

// Midpoint-radius enclosure: when status is ok or unbounded, true value in [mid - rad, mid + rad].
template <class REAL>
struct Rmidrad {
    REAL    mid;     // approximation (finite when status == ok)
    REAL    rad;     // >= 0; +inf when status == unbounded
    Rstatus status;
};

} // namespace vmplapack
```

Note: `non_finite` is **only** for non-finite *inputs*. A non-finite value *produced* from finite
inputs is overflow → `unbounded` (see §6.3). Do not re-add "or NaN produced" to this comment.

---

## 5. Public API (`include/vmplapack/vmpblas.h`) ― BLAS-consistent argument order

```cpp
namespace vmplapack {

// ---- accurate (compensated; bare REAL). Boundary/empty rules in §8. ----
template <class REAL>
REAL Rsum(std::ptrdiff_t n, const REAL* x, std::ptrdiff_t incx = 1);

template <class REAL>
REAL Rdot(std::ptrdiff_t n,
          const REAL* x, std::ptrdiff_t incx,
          const REAL* y, std::ptrdiff_t incy);

// ---- verified (rigorous midrad enclosure). Boundary/empty rules in §8. ----
template <class REAL>
Rmidrad<REAL> vRsum(std::ptrdiff_t n, const REAL* x, std::ptrdiff_t incx = 1);

template <class REAL>
Rmidrad<REAL> vRdot(std::ptrdiff_t n,
                    const REAL* x, std::ptrdiff_t incx,
                    const REAL* y, std::ptrdiff_t incy);

// ---- verified residual r = b - A x ; row-major general A (m x n), lda >= n ----
// Writes one Rmidrad per row into out[0..m-1]; returns the worst Rstatus (Rstatus_rank order).
template <class REAL>
Rstatus vRresidual(std::ptrdiff_t m, std::ptrdiff_t n,
                   const REAL* A, std::ptrdiff_t lda,
                   const REAL* x,
                   const REAL* b,
                   Rmidrad<REAL>* out);

} // namespace vmplapack
```

Stride convention (BLAS): for `incx >= 1`, element `i` is `x[i*incx]`. The accurate/verified API
asymmetry (bare `REAL` vs `Rmidrad`+`Rstatus`) is deliberate; only verified routines carry a status
channel. All empty/invalid/boundary behavior is defined in §8.

---

## 6. Preconditions and guarantees (read before §7)

### 6.1 EFT operand-aliasing safety (CRITICAL)

The reduction pseudocode reuses an accumulator as both input and output, e.g.
`Rtwosum(s, x[i]) -> (s, e)` and `Rtwosum(p, h) -> (p, q)`. The natural C++ calls

```cpp
Rtwosum(s, x[ix], s, e);     // output s aliases input
Rtwosum(p, h,     p, q);     // output p aliases input
```

**must be correct.** Therefore EFTs take their input operands **by value** (never `const REAL&` read
after an aliased output is written). The MPFR copy cost is accepted ("correctness over cleverness").

### 6.2 EFT exactness preconditions

`Rtwosum` reconstructs `a+b = s+e` *exactly* when inputs are finite, the sum does not overflow,
subnormals are not flushed (§9), and rounding is round-to-nearest. `Rfasttwosum` additionally
requires `abs(a) >= abs(b)`.

`Rtwoproduct` reconstructs `a*b = p+e` *exactly* under the same conditions **and** when the exact
residual `a*b - fl(a*b)` is representable in the working format ― i.e. when the leading product does
**not** overflow **and does not underflow**. Under product underflow, `e = fma(a,b,-p)` is the
*correctly-rounded* residual but not necessarily the *exact* one. Therefore: tests that assert exact
`hi+lo` reconstruction must avoid product underflow as well as overflow. Underflow / subnormal cases
are tested as **robustness** cases (correct status, valid enclosure), **not** as exact-TwoProduct
cases, unless a theorem + implementation path covering underflow is explicitly added.

### 6.3 Verified certificate semantics (NaN / Inf classification ― single source of truth)

Verified routines must **first scan all inputs**. If any input is NaN or Inf, return
`status == non_finite` and **claim no certificate** (the true value is not a finite real). Use
`is_finite` for this scan (it is false for both NaN and Inf ― no NaN/Inf distinction is needed).

After all inputs are confirmed finite, the true real value is finite, so **any** non-finite endpoint,
midpoint, or radius produced by the directed passes or by midrad construction is overflow / loss of
finite bounds: return `status == unbounded` with `mid = 0`, `rad = +inf` ― the valid but useless
enclosure `[-inf, +inf]`. **Never claim finite verification for an overflowed native computation.**

`Rmake_midrad` (§7.3) assumes input non-finiteness has **already** been classified by the caller; it
is **not** an input classifier. It only turns already-classified directed bounds into a midrad
enclosure and maps any residual non-finiteness/invariant breach to `unbounded`.

---

## 7. Algorithms

### 7.1 EFT (`vmplapack_eft.h`) ― by value, default round-to-nearest

```cpp
// a + b = s + e exactly,  s = fl(a + b)   (Knuth TwoSum, branch-free, alias-safe)
template <class REAL>
void Rtwosum(REAL a, REAL b, REAL& s, REAL& e) {
    s = a + b;
    REAL z  = s - a;
    REAL t1 = s - z;
    REAL t2 = a - t1;
    REAL t3 = b - z;
    e = t2 + t3;
}

// Precondition: abs(a) >= abs(b), and a + b does not overflow.
template <class REAL>
void Rfasttwosum(REAL a, REAL b, REAL& s, REAL& e) {
    s = a + b;
    REAL z = s - a;
    e = b - z;
}

// a * b = p + e exactly (under §6.2),  p = fl(a * b)   (FMA-based; alias-safe)
template <class REAL>
void Rtwoproduct(REAL a, REAL b, REAL& p, REAL& e) {
    p = a * b;
    REAL neg_p = -p;
    e = Rarith<REAL>::fma(a, b, neg_p);    // std::fma / fmaf / mpfrxx::fma
}
```

Provide a Dekker `Rsplit` + split-based `Rtwoproduct` under `VMPLAPACK_USE_DEKKER_TWOPRODUCT`. Split
factor `2^ceil(p/2)+1`: float `0x1p12f+1`, double `0x1p27+1`, mpfr `2^ceil(W/2)+1`. **Exact only when
`factor*a` does not overflow and subnormals are not flushed**; extreme-magnitude tests must avoid
split overflow or expect `unbounded`/`non_finite`.

### 7.2 `Rsum` ― Sum2; `Rdot` ― Dot2 (OgitaRumpOishi), guarded + fully materialized

```text
Rsum:
    if n <= 0: return zero()
    s = x[0]; c = zero()
    for i = 1..n-1:  Rtwosum(s, x[i], s, e);  c = c + e
    REAL res = s + c;  return res
```

```text
Rdot:
    if n <= 0: return zero()
    Rtwoproduct(x[0], y[0], p, s)
    for i = 1..n-1:
        Rtwoproduct(x[i], y[i], h, r)
        Rtwosum(p, h, p, q)
        REAL t  = q + r          // expand: one op per statement (materialization rule)
        REAL sn = s + t
        s = sn
    REAL res = p + s
    return res
```

### 7.3 Midrad helpers (`vmplapack_utils.h`) ― shared by M5 and M6

```cpp
#include <cassert>

// Overflow-safe midpoint of [lo, hi]: lo/2 + hi/2 never overflows for finite lo, hi.
template <class REAL>
REAL Rmidpoint(REAL lo, REAL hi) {
    REAL h = Rarith<REAL>::half();
    REAL a = lo * h;
    REAL b = hi * h;
    REAL m = a + b;
    return m;
}

// |a - b| rounded UPWARD (over-estimates the true distance).
template <class REAL>
REAL Rupward_abs_diff(REAL a, REAL b) {
    typename Rarith<REAL>::round_up g;
    REAL d = (a >= b) ? (a - b) : (b - a);
    return d;
}

// Convert already-classified directed bounds lo <= true <= hi (finite INPUTS guaranteed by caller)
// plus a chosen midpoint `mid` into a rigorous midrad enclosure. NOT an input classifier (§6.3).
template <class REAL>
Rmidrad<REAL> Rmake_midrad(REAL lo, REAL hi, REAL mid) {
    using A = Rarith<REAL>;
    if (!A::is_finite(lo) || !A::is_finite(hi)) {            // a bound overflowed
        return { A::zero(), A::infinity(), Rstatus::unbounded };
    }
    if (hi < lo) {                                            // broken invariant (contract breach)
        assert(false && "Rmake_midrad: lo > hi");
        return { A::zero(), A::infinity(), Rstatus::unbounded };
    }
    REAL m = mid;
    if (!A::is_finite(m)) {                                   // bounds finite but midpoint failed:
        m = Rmidpoint(lo, hi);                               // keep the finite enclosure
    }
    REAL r1  = Rupward_abs_diff(m, lo);
    REAL r2  = Rupward_abs_diff(hi, m);
    REAL rad = (r1 > r2) ? r1 : r2;
    if (!A::is_finite(rad)) {
        return { A::zero(), A::infinity(), Rstatus::unbounded };
    }
    return { m, rad, Rstatus::ok };
}
```

`Rmake_midrad` is rigorous even when `mid` lies outside `[lo, hi]` (one upward distance covers the
whole interval), and never discards a finite directed enclosure merely because the accurate midpoint
failed.

### 7.4 `vRdot` ― verified via directed-rounding enclosure

```text
1. Validate per §8 (n<0/inc/pointer -> invalid_input; n==0 -> {zero,zero,ok}).
2. If any x[i] or y[i] is non-finite -> return { 0, +inf, non_finite }.   // no certificate (§6.3)
3. Upper bound (rounded toward +inf):
       { round_up   g;  acc = zero(); for i: { prod = x[ix]*y[iy]; acc = acc + prod; }  sup = acc; }
4. Lower bound (rounded toward -inf):
       { round_down g;  acc = zero(); for i: { prod = x[ix]*y[iy]; acc = acc + prod; }  inf = acc; }
   // true dot in [inf, sup]
5. mid = Rmidpoint(inf, sup)                 // ambient round-to-nearest
6. return Rmake_midrad(inf, sup, mid)        // ok / unbounded handled inside
```

Plain product + add only (no FMA) so both passes are structurally identical on every tier.
Materialize each `prod` and each `acc`. `vRsum` is the same without the product step.

### 7.5 M6 upgrade ― accurate *and* verified

Same `[inf, sup]` enclosure; sharpen the midpoint with Dot2:

```text
mid = Rdot(x, y)                 // accurate midpoint, ambient round-to-nearest
return Rmake_midrad(inf, sup, mid)   // if mid is non-finite while [inf,sup] finite, the helper
                                     // falls back to Rmidpoint(inf,sup) ― finite enclosure kept
```

Still rigorous (§7.3), now accurate, with no analytical constant.

*(Optional)* Also implement the ORO a-priori Dot2 bound `|res － xy|  u|res| + γΣ|xy|` and
compare tightness. Two hard requirements: take the **exact** `γ` (subscript + validity `nu < 1`)
from ORO 2005; and `Σ|xy|` must be a rigorous **upper bound** ― compute it by an **upward-directed
pass** (or take the upper endpoint of `vRdot(|x|, |y|)`), **never** via accurate `Rdot`. Round the
whole expression outward.

### 7.6 `vRresidual` ― r = b - A x (sign-correct)

```text
Validate per §8 (m<0/n<0 -> invalid_input; m==0 -> ok, write nothing; m>0,n==0 -> r_i = b_i).
For each row i:
  if b_i or any A[i,j] or x[j] (0<=j<n) is non-finite: out[i] = {0,+inf,non_finite}; continue
  upper (RNDU):  acc = b_i;  for j: { na = -A[i,j];  prod = na*x[j];  acc = acc + prod; }  sup = acc
  lower (RNDD):  acc = b_i;  for j: { na = -A[i,j];  prod = na*x[j];  acc = acc + prod; }  inf = acc
  mid = Rmidpoint(inf, sup);  out[i] = Rmake_midrad(inf, sup, mid)
Return the worst status by Rstatus_rank over all written rows.
```

Each pass starts from `b_i` and adds `(-A[i,j])*x[j]` (NOT `-b_i` plus `Ax`). Example-level only;
do **not** build a solver.

---

## 8. Boundary conditions and input validation

```text
Rsum / Rdot (accurate, bare REAL):
- n <= 0 returns zero() and inspects no pointer.
- If n > 0: incx, incy >= 1 and all accessed pointers valid are PRECONDITIONS.
  Violation is UB in release and an assertion failure in debug builds.

vRsum / vRdot (verified):
- n < 0  -> { 0, +inf, invalid_input }.
- n == 0 -> { zero, zero, ok } and inspects no pointer.
- n > 0  -> incx, incy >= 1 and non-null accessed pointers required; otherwise invalid_input.

vRresidual (verified):
- m < 0 or n < 0 -> invalid_input (write nothing).
- m == 0 -> ok, write nothing; pointers may be null.
- m > 0, n == 0  -> residual is r_i = b_i for each row; b and out must be valid; A and x are not read.
- m > 0, n > 0   -> A, x, b, out valid and lda >= n required; otherwise invalid_input.
```

These rules are normative; tests must cover each case.

---

## 9. Rounding / floating-point contract (mandatory)

### 9.1 Native tiers (`float`, `double`)

For **every** TU with EFT / verified / native `Rarith` code (centralize in `RoundingControl.cmake`,
PUBLIC/INTERFACE so downstream TUs inherit them):

```text
-ffp-contract=off     # no implicit a*b+c -> fma
-fno-fast-math
-frounding-math       # the verified path changes the rounding mode at runtime
```

(MSVC `/fp:strict`.) `-frounding-math` is documented by GCC as experimental and does not by itself
guarantee disabling every rounding-affected optimization, so **also** add, in each such TU:

```cpp
#pragma STDC FENV_ACCESS ON
#pragma STDC FP_CONTRACT OFF
```

and put this guard in the umbrella header (so any fast-math downstream build fails to compile ―
this compile-time rejection is the PRIMARY required safeguard):

```cpp
#if defined(__FAST_MATH__)
#error "vMPLAPACK requires IEEE-754 semantics; do not compile with -ffast-math."
#endif
```

Require `FLT_EVAL_METHOD == 0` (no x87 excess precision ― also what makes the `float` tier sound,
since float ops must round to binary32); on 32-bit x86 add `-mfpmath=sse -msse2`:

```cpp
#include <cfloat>
static_assert(FLT_EVAL_METHOD == 0, "Excess precision breaks EFT, interval bounds, and the float tier.");
```

**Subnormals must not be flushed.** `-ffast-math` (and sometimes other settings) enable FTZ/DAZ in
the SSE MXCSR, which breaks EFT and the enclosure. Ensure FTZ/DAZ are OFF and verify it at runtime
(§9.3c).

### 9.2 MPFR tier

MPFR is correctly rounded by construction; directed rounding via `mpfrxx::rounding_mode_scope`. Do
**not** rely on MPFR FMA auto-fusion: `Rtwoproduct` calls `mpfrxx::fma` explicitly; the verified path
uses plain product + add. The accurate path runs in `MPFR_RNDN`.

### 9.3 Environment self-tests (`test_Rarith_environment.cpp`), per tier

a) the rounding mode returns to nearest after a `round_up`/`round_down` scope exits and on exception;
b) crafted cancellations reconstruct exactly via `Rtwosum`/`Rtwoproduct` (under §6.2 preconditions);
c) **subnormals are alive** (the smallest positive subnormal is not flushed) **and** entering
   `round_up`/`round_down` **actually changes** the result of an ordinary `a + b` and `a * b` at
   runtime (proving the directed-rounding contract is in force, not optimized away);
d) the directed passes bracket a known hard case.

Negative fast-math test: the PRIMARY required behavior is **compile-time rejection** via the
`__FAST_MATH__` guard. Optionally also add a hand-crafted inclusion test that fails **if the guard is
intentionally disabled**. Do **not** require arbitrary `-ffast-math` builds to fail numerically ― they
may not, and not failing does not imply soundness.

---

## 10. MPFR fixed precision `W` discipline + oracle (`tests/Rdot_oracle.h`, TEST ONLY)

### 10.1 Fixed `W` discipline

MPFR's default precision affects only **subsequently** initialized variables; previously created
objects keep their precision. Therefore:

```text
- W is fixed per test process. Call mpfrxx::set_default_precision_bits(W) BEFORE constructing any
  mpfr_class used as kernel input. Changing W after any mpfr_class exists is forbidden.
- Every mpfr_class input must have precision W. Add a test helper asserting mpfr_get_prec(x.mpfr_data()) == W.
- Run W=53 and W=512 as SEPARATE CTest cases (e.g. MPFRXX_DEFAULT_PRECISION_BITS=53 / =512), or,
  within one executable, set W and fully reconstruct all MPFR objects before each block.
```

### 10.2 Oracle = naive computation at precision `P`, returned as a directed INTERVAL

Compute the reductions **naively** (no EFT) at a precision `P` far above the tier under test, and
return a rigorous interval, not a point:

```text
ref_lo = naive dot at precision P, MPFR_RNDD     // <= true
ref_hi = naive dot at precision P, MPFR_RNDU     // >= true
```

`P` selection: start from `P = max(512, 4*precision_bits(tier) + 64)` as a heuristic; the **rigorous
safeguard is precision doubling** ― also compute `[ref_lo, ref_hi]` at `2P` and require the two
intervals to agree to far below the tier's `u`. Large exponent spans / deep cancellation can need
more than the heuristic; widen `P` until the doubling check passes. Materialize each oracle op into
`mpfr_class`; keep MPFR FMA auto-fusion off.

Exact input widening into the oracle (oracle must see the *same* numbers):

- `REAL = float`/`double` → `mpfrxx::mpfr_class::with_precision(P, x[i])` (exact for `P >= p`).
- `REAL = mpfr_class@W` → `mpfrxx::mpfr_class(x[i].mpfr_data(), P)` (exact for `P >= W`). Do **not**
  move across differing precisions; construct/copy at `P`.

---

## 11. Tests

Lightweight framework of your choice. Every test states the failure it catches. Run the suite for
**all four** instantiations: `float`; `double`; `mpfr_class @ W=53` (mirrors binary64 precision ―
cross-check at equal precision; unbounded exponent range, NOT bit-exact binary64); `mpfr_class @
W=512` (primary working precision).

### 11.1 Accurate kernels ― condition-aware (do NOT use `|got-ref| <= C*u*|ref|`)

A bare relative test breaks under cancellation / exact zero. Use the Dot2 error shape with a rigorous
**upper bound** on the term sum:

```text
S_hi   = oracle upper bound of Σ|x_i y_i| via an UPWARD-directed high-precision pass
ref_lo, ref_hi = oracle interval (§10.2);  ref = midpoint of [ref_lo, ref_hi]
require |got - ref| <= u_tier*|ref| + C * u_tier^2 * S_hi + floor
        // C ~ 4 ; floor = a few * (ref_hi - ref_lo)  (the oracle interval width; tier-agnostic)
```

This passes Dot2 yet fails a non-compensated dot (error ~ `u*S` ≫ `u^2*S` when `S ≫ |ref|`). The
criterion is an explicit **test heuristic** (constant `C`) unless the exact ORO theorem constant is
implemented; if the oracle interval width is not negligible vs the right-hand side, increase `P`. Do
**not** use a native "smallest-normal" floor for the mpfr tier ― the oracle-width floor above is
correct for every tier.

### 11.2 Verified inclusion ― enclose the oracle INTERVAL (core tests)

```text
m = vRdot(...)
require m.status == ok or unbounded
require m.rad >= 0
require (m.mid - m.rad) <= ref_lo  AND  ref_hi <= (m.mid + m.rad)     // in mpfr_class, exact widening
```

Must fail if the verified interval does not cover `[ref_lo, ref_hi]`; `rad < 0`; `status == ok` with a
non-finite produced; or a fast-math/contraction build changes the result. For deliberately overflowing
inputs require `status == unbounded` (valid `[-inf,+inf]`); for non-finite inputs require
`non_finite` ― never a false finite certificate.

### 11.3 Case taxonomy (each kernel, each tier)

Empty (n==0), singleton, all-positive, alternating signs, large exponent gaps, tiny/subnormal
magnitudes (with subnormals confirmed alive, §9.3c), exact cancellation to zero, overflow-inducing
magnitudes (expect `unbounded`), non-finite inputs (expect `non_finite`), invalid n/stride/pointer
(expect `invalid_input` for verified; debug-assert for accurate), random (fixed seed). Strides
`incx==1` and non-unit (e.g. 3).

### 11.4 Ill-conditioned generators (`tests/Rgendot.h`)

M8 closes the earlier random-generation deferral. `Rgendot.h` provides both deterministic adversarial
families and a seeded randomized high-condition generator. All generated input values must be exactly
constructed in the target tier: native `float`/`double` use in-range `std::ldexp`; mpfr uses
`mpfr_set_ui_2exp(_, 1, j, MPFR_RNDN)` at precision `W`; no generator constructs tier inputs through
`double`.

Deterministic families:

```text
Family A (alternating):  x_i = 1; y chosen so the exact dot = delta.
Family B (two-term):     x = [1, 1], y = [1, y2], y2 exactly representable; exact dot = 1 + y2.
Family C (exponent-cancellation):
    x = [2^j, s, 2^j], y = [1, 1, -1]
    The two 2^j terms cancel exactly; placing s between them makes naive sequential accumulation lose s
    when j >= W. true dot = s.
M8 adversarial families:
    heavy cancellation, alternating signs, huge/small scale mixing, and exact nontrivial cancellation.
Ordering variants:
    generated, sorted, reversed, shuffled variants of the same multiset.
```

Condition-number convention is fixed and both values are recorded:

```text
cond_oro = 2 * S / |dot|        # ORO 2005 convention; use for target/benchmark axes
cond_sum =     S / |dot|        # sum-of-absolute-products convention
S = sum_i |x_i*y_i|
```

For a nonzero exact dot, `cond_sum >= 1` and `cond_oro >= 2`. A finite `target_cond < 2` is invalid
under the ORO convention. If exact `dot == 0` and `S > 0`, record
`cond_oro = cond_sum = +inf`; that is a valid adversarial case but is excluded from finite log-scale
target matching.

Each generated case records:

```text
target_cond, measured_cond_mpfr (= cond_oro from the oracle), measured cond_sum,
seed, scale, permutation, tier, status
```

The measured condition uses the §10.2 MPFR oracle interval and an upward high-precision pass for `S`,
not a bare W-bit mpfr evaluation. A targeted finite generator must either satisfy:

```text
|log10(measured_cond_mpfr) - log10(target_cond)| <= 0.25 decades
```

or report `unachievable` for that tier, length, exponent range, and target. Tests should use several
finite target conditions per tier, spanning low, medium, and high condition regimes, so oracle comparison
is not tied to a single hand-picked condition number.

Generator tests must check the theoretical scaling variables, not a crude `cond ~ 1/u` threshold:

```text
naive dot scale:  about gamma_n * cond_sum, equivalently 0.5 * gamma_n * cond_oro
Dot2/Rdot scale: about u + 0.5 * gamma_n^2 * cond_oro
```

Across seeds, multiple target conditions, deterministic families, and ordering variants, `vRdot` must
always enclose the oracle interval. For every finite targeted high-condition random case,
`vRdot_apriori` must also be compared with the oracle interval and enclose it with `status == ok`.
`Rdot` may be inaccurate for high conditions; verified inclusion must not fail.

### 11.5 Cross-tier consistency

Restrict to inputs and intermediate results well within normal binary64 range, no subnormals, no
overflow. Under those restrictions, `Rdot<double>` and `Rdot<mpfr_class@W=53>` must agree exactly or
within a small **result-ulp** tolerance (ulp of the result magnitude, not of `2^-53`), and their
`vRdot` enclosures must overlap. Outside those restrictions, require only that the verified enclosures
overlap. (The exponent-range / subnormal differences between binary64 and unbounded-range mpfr are
why bit-exact equality is not required; clamping mpfr's exponent range is out of scope.)

### 11.6 EFT exactness tests ― via MPFR high precision, not same-tier `hi + lo`

Computing `hi + lo` in the same tier re-rounds it. For random exact-reconstruction tests, widen
`a, b, hi, lo` exactly to oracle precision `P` and check in MPFR that `exact(ab) == exact(hi) +
exact(lo)` for ` ∈ {+, *}`. Crafted cases with a known decomposition may instead assert the exact
`hi`/`lo` constants directly (no MPFR). In native-only builds (MPFR off), the MPFR-based exact tests
are skipped; the crafted-constant tests still run.

---

## 12. CMake

Options:

```text
VMPLAPACK_ENABLE_MPFR                 (default ON)    # mpfr tier in the library AND MPFR oracle in tests
VMPLAPACK_ENABLE_TESTS               (default ON)
VMPLAPACK_ENABLE_EXAMPLES            (default ON)
VMPLAPACK_USE_DEKKER_TWOPRODUCT      (default OFF)
VMPLAPACK_ENABLE_SANITIZERS          (default OFF)   # ASan+UBSan instrumentation for Debug validation
VMPLAPACK_ENABLE_FAST_MATH_GUARD_TEST (default ON)    # try_compile must reject -ffast-math via __FAST_MATH__
```

- `cmake/RoundingControl.cmake` defines an INTERFACE target carrying the §9.1 flags; link it into the
  library, tests, examples, and propagate it **PUBLIC/INTERFACE** so consumers inherit the IEEE-754
  requirements. README must instruct users to link `vmplapack::vmplapack`, not include headers with
  arbitrary flags.
- `VMPLAPACK_ENABLE_MPFR=ON` → find GMP + MPFR (+ MPC if using combined `gmpfrxx_mkII.h`) and add
  gmpfrxx_mkII's `include/` to the library + tests; mpfr tier + oracle compile. `OFF` → native-only
  (float+double), **no** GMP/MPFR anywhere; oracle-based and MPFR-exactness tests are **skipped**, and
  only the environment + crafted-constant native-EFT smoke tests are built (state this explicitly).
- C++17; `-Wall -Wextra -Werror` recommended, never relaxing the FP flags. `enable_testing()` +
  `add_test`; register the W=53 and W=512 mpfr suites as separate tests (§10.1).
- M9 configuration must run the fast-math negative compile check by default: a `try_compile` using
  `-ffast-math` includes the umbrella header and must fail through the `__FAST_MATH__` `#error`.
- `VMPLAPACK_ENABLE_SANITIZERS=ON` adds ASan+UBSan instrumentation for Debug validation; it must not
  relax the strict FP flags.

`cmake -S . -B build && cmake --build build && ctest` must pass with GMP/MPFR present.

---

## 13. README.md

Cover: vMPLAPACK is a lightweight MPLAPACK-like prototype, **not** part of MPLAPACK; the MPLAPACK
idiom; the three co-equal tiers and that mpfr is first-class with mpfr-at-high-precision additionally
serving as the oracle; accurate vs verified and the inclusion contract, including the
`unbounded`/`non_finite` certificate semantics (§6.3); the boundary rules (§8); the rounding/FP
contract (§9) and **why** it is mandatory (consumers must link the target); the fixed-`W` discipline;
how to build/run `ctest`. Reading list with identifiers:

```text
- T. Ogita, S. M. Rump, S. Oishi, "Accurate Sum and Dot Product",
  SIAM J. Sci. Comput. 26(6), 19551988, 2005. DOI: 10.1137/030601818.
- S. M. Rump, "Verification methods: Rigorous results using floating-point arithmetic",
  Acta Numerica 19, 287449, 2010. DOI: 10.1017/S096249291000005X.
- N. J. Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., SIAM, 2002.
  ISBN: 0-89871-521-0 / 978-0-89871-521-7. DOI: 10.1137/1.9780898718027.
```

---

## 14. Milestones (one PR each; pass acceptance before advancing)

- **M0 ― Skeleton + native arith.** Layout; CMake + `RoundingControl.cmake` (PUBLIC flags, pragmas,
  `__FAST_MATH__` guard); `vmplapack_utils.h` (`Rstatus` incl. `unbounded`; `Rstatus_rank`; `Rmidrad`;
  `Rmidpoint`, `Rupward_abs_diff`, `Rmake_midrad` with lo>hi defensive check and midpoint fallback);
  `Rarith<float>`/`Rarith<double>` (signed `precision_bits`, `half/one/infinity`, non-copyable/movable
  scopes); `test_Rarith_environment.cpp` (§9.3 ad). *Accept:* builds; environment test green.
- **M1 ― EFT, native.** By-value `Rtwosum`/`Rfasttwosum`/`Rtwoproduct` (+ Dekker `Rsplit`).
  Exact-reconstruction tests **including aliased calls** `Rtwosum(s,x,s,e)`, `Rtwosum(p,h,p,q)`, under
  §6.2 preconditions, for float and double, both TwoProduct paths; crafted-constant exact checks run
  native-only, random exact checks via MPFR widening (§11.6) when MPFR is enabled. *Accept:* exact
  `ab=hi+lo` under §6.2, alias-safe.
- **M2 ― MPFR tier + API probe.** First a compile/configure probe verifying `mpfrxx::mpfr_class`,
  `mpfrxx::fma`, `mpfrxx::rounding_mode_scope` with `MPFR_RNDU/RNDD`, `mpfr_data()`/`get_mpfr_t()`, and
  a precision-setting constructor exist; **stop and report** the exact missing symbol if not. Then
  `Rarith<mpfrxx::mpfr_class>` (signed-cast 2^-W); re-run M1 EFT tests for mpfr at W=53 and W=512;
  assert `mpfr_get_prec == W`. *Accept:* probe passes; EFT exact on mpfr; all tiers co-equal.
- **M3 ― Oracle + generators.** `Rdot_oracle.h` returning `[ref_lo,ref_hi]` with precision-doubling
  self-check and exact widening; deterministic `Rgendot.h` (A/B/C, powers exact in-tier). *Accept:*
  oracle interval stable under `P→2P`; widening exact.
- **M4 ― Accurate.** `Rsum`/`Rdot` (guarded, materialized) + explicit instantiations + `extern
  template`. Condition-aware tests (§11.1) on all tiers; cross-tier consistency (§11.5). *Accept:*
  condition-aware bound holds; double vs mpfr@53 agree under the restricted regime.
- **M5 ― Verified + residual.** `vRsum`/`vRdot` (directed enclosure via `Rmake_midrad`; non-finite
  inputs→`non_finite`; overflow→`unbounded`; §8 validation); sign-correct `vRresidual` returning the
  worst `Rstatus`; examples. Inclusion-of-`[ref_lo,ref_hi]` tests across the sweep incl. Family C
  `cond > 1/u`, plus explicit overflow, non-finite, and invalid-input cases, on all tiers. *Accept:*
  inclusion never fails; statuses correct.
- **M6 ― Accurate + verified.** `mid = Rdot(...)` via `Rmake_midrad` (finite-enclosure fallback when
  the accurate midpoint is non-finite). *(Optional)* ORO a-priori bound with `Σ|x_iy_i|` as a rigorous
  **upper** bound (upward pass), tightness comparison. README finalized. *Accept:* midpoint accuracy
  improves; inclusion still always holds.

---

## 15. Global acceptance criteria

```text
1.  Builds independently with CMake; no real-MPLAPACK dependency.
2.  Native-only build (float+double) needs no GMP/MPFR; ENABLE_MPFR=ON builds the mpfr tier and oracle.
3.  EFTs are alias-safe; exact reconstruction holds under §6.2 (no product overflow/underflow), verified
    via MPFR high precision (not same-tier hi+lo).
4.  Rsum/Rdot are accurate per the condition-aware criterion (S as an upward upper bound) on ALL tiers;
    relative-only tests are NOT used.
5.  vRdot encloses the oracle interval [ref_lo, ref_hi] for ALL tested inputs on ALL tiers, incl. Family C
    with cond > 1/u_tier.
6.  NaN/Inf INPUT -> non_finite (no certificate); finite-input overflow -> unbounded with rad=+inf (valid
    [-inf,+inf]); no false finite certificate is ever produced. NaN/Inf classification has a single source
    of truth (§6.3); Rmake_midrad is not an input classifier.
7.  Empty/invalid/boundary behavior matches §8 exactly (n==0, n<0, null pointers, vRresidual m==0 / n==0).
8.  vRresidual computes r = b - A x (from b_i, adding (-A)*x) and returns the worst Rstatus by Rstatus_rank.
9.  The rounding/FP contract is enforced per tier and verified at runtime (mode restore, subnormals alive,
    directed rounding affects +/*), exposed PUBLIC; __FAST_MATH__ is rejected at compile time.
10. MPFR W is fixed; inputs constructed after setting W and verified to have precision W; 2^-W uses a signed
    exponent; W=53 and W=512 run as separate test cases.
11. Public names follow MPLAPACK style with the scalar type as the REAL template parameter; accurate
    routines return REAL, verified routines carry Rstatus.
```

---

## 16. Agent operating rules

- Implement strictly in milestone order; one PR per milestone; pass acceptance before advancing.
- Do not alter the public signatures, `Rmidrad`/`Rstatus`, the `Rarith<REAL>` interface, the
  preconditions/guarantees (§6), the boundary rules (§8), or the rounding contract (§9) without
  stopping to report.
- Keep `float`, `double`, `mpfrxx::mpfr_class` strictly co-equal: one generic template per kernel; all
  tier differences live in `Rarith<REAL>`.
- Soundness over tightness, correctness over cleverness. All code, identifiers, and comments in English.

---
# docs/SPEC.md §14 — addendum: milestones M7–M11

Paste into `docs/SPEC.md §14` after M6. These consolidate the M0–M6 prototype into a shippable kernel
library. **All of M7–M11 remain standalone (no MPLAPACK dependency).** Same rules as M0–M6: one PR per
milestone; pass acceptance before advancing; do not change frozen contracts (§5/§6/§8/§9) without
stopping to report.

---

- **M7 — ORO a-priori accuracy bound.** Promote the M6-optional bound to a routine. Implement
  `vRdot_apriori<REAL>` returning an `Rmidrad<REAL>` whose `mid` is the Dot2 result (`Rdot`) and whose
  `rad` is a **rigorous, outward-rounded, computable** bound derived from the Ogita–Rump–Oishi Dot2
  a-priori estimate.

  Definitions (fix precisely):
  - `u` = **unit roundoff** = `2^-p` = ½ ulp(1), `p = Rarith<REAL>::precision_bits()`. (NOT machine
    epsilon `2^-(p-1)`.)
  - `γ_n = n*u / (1 − n*u)`, with the ORO validity condition `n*u < 1` (n = vector length).
  - `S = Σ |x_i*y_i|`; compute `S_up >= S` with the **named primitive**
    `Rsum_abs_dot_upper<REAL>(n, x, incx, y, incy) -> REAL` — in an **upward-rounding** scope, accumulate
    `|x[i]|*|y[i]|` with **upward-rounded products and summation** (native `FE_UPWARD` / mpfr
    `MPFR_RNDU`). It **must not** call accurate `Rdot`. It returns a bare `REAL` and yields `+inf` on
    overflow; finite inputs are pre-checked by the caller (§6.3), so `+inf` here means overflow and is
    mapped to `unbounded`.

  The bound. The raw ORO Dot2 estimate is `|res − xᵀy| ≤ u|xᵀy| + γ_n²·S` (no underflow). Since `res`,
  not `xᵀy`, is what we have, substitute `|xᵀy| ≤ |res| + |e|` and solve, giving the **computable**
  form. A valid implementation:

  ```text
  rad >= outward( ( u*abs(mid) + γ_n^2 * S_up + underflow_term ) / (1 - u) )   # every op rounded outward
  ```

  - `underflow_term`: the **ORO 2005 Prop. 5.5 underflow term, `5*n*η`** (η = the underflow unit). For
    native float/double `η = Rarith::eta() = denorm_min`. **For the mpfr tier, do NOT claim that the
    wide exponent range prevents underflow** — products of very small values can fall below MPFR's
    exponent lower bound. The implementation must either **(a)** check the MPFR underflow flag /
    exponent condition for the bound pass and return `unbounded` (or include the term) if underflow
    occurred, or **(b)** enforce and **test** a no-underflow precondition. (M7 adds `Rarith<REAL>::eta()`
    — a sanctioned extension of §2; for mpfr it is the smallest representable, but the bound path relies
    on (a)/(b), not on `eta` alone.)
  - `1 − u > 0` always; round the division outward.

  Gate. Evaluate `n*u < 1` **conservatively with integer arithmetic, never in floating point near the
  boundary**. Since `u = 2^-p`, the exact condition is `n < 2^p`. Let `length_type` be the length
  parameter's type (e.g. `mplapackint`/`std::ptrdiff_t`; a signed type's `digits` excludes the sign
  bit). If `p >= std::numeric_limits<length_type>::digits`, then `2^p` exceeds every representable
  length and the gate **passes** for all `n` — and `2^p` is not representable in the type, so do not try
  to form it (this is the correct boundary: at `p == digits`, `length_type{1} << p` would be UB).
  Otherwise compare `n < (length_type{1} << p)` exactly. If the gate cannot be decided exactly, fail
  conservatively. When it fails (`n*u ≥ 1`), return `{ Rdot(x,y), +inf, unbounded }` — keep the Dot2
  midpoint (still a useful approximation); only the radius is meaningless.

  Disambiguation (write this in the code/SPEC): **this is the a-priori Dot2 estimate, NOT ORO Algorithm
  5.8 `Dot2Err`.** `Dot2Err` is a separate pure-floating-point error-bound routine with a different
  computable bound and a stronger `2*n*u < 1` gate. Do not implement one as the other.

  Non-finite policy (matches the verified kernels, §6.3/§8): non-finite **input** → `non_finite`;
  finite-input overflow in `res`, `S_up`, or the bound → `unbounded`; signed zeros via IEEE
  (`Σ|·|` maps `−0`→`+0`).

  Tie-in: once M7 lands, the M4 accurate-kernel test (§11.1) may be **strengthened** from the heuristic
  `|got−ref| ≤ u|ref| + C·u²·S` to the rigorous `|got − ref| ≤ vRdot_apriori.rad` (the heuristic's fixed
  `C` does not absorb the `γ_n²` n-growth; the M7 bound does).

  *Tasks:* `Rsum_abs_dot_upper` (with the `+inf`-on-overflow contract); the outward bound with the
  `/(1−u)` correction; the integer `n*u<1` gate; explicit underflow handling (`5*n*η` term **or** checked
  no-underflow precondition) with `Rarith::eta()`; tightness comparison vs `vRdot.rad`. *Accept:*
  `rad ≥ |mid − ref|` for ALL tested inputs on ALL tiers (incl. cancellation and, separately, underflow
  cases); `Rsum_abs_dot_upper` proven `≥ S` against the oracle; the integer gate behaves correctly at
  the boundary; tightness vs `vRdot.rad` reported. **Must not** use accurate `Rdot` for `S_up`.

- **M8 — Random + adversarial high-condition generators.** Close the deferred random `gendot` (§11.4)
  **and** add deterministic adversarial families. Implement a **seeded** Ogita–Rump–Oishi-style
  randomized generator targeting a prescribed dot condition number, **and** deterministic adversarial
  generators: heavy cancellation, alternating signs, huge/small scale mixing, and **input ordering
  variants** (sorted, reversed, shuffled) of the same multiset. All tier values constructed **exactly
  in-tier** per the **§10.2 exact-construction rule** (mpfr via `mpfr_set_ui_2exp`, native via in-range
  `ldexp`; never through `double`).

  Condition-number convention (fix it — the two differ by a factor 2): record **both**

  ```text
  cond_oro = 2 * S / |dot|        # the ORO 2005 definition (use this for acceptance/benchmark axes)
  cond_sum =     S / |dot|        # the "sum of |x_i y_i|" convention
  ```

  Because `S = Σ|x_iy_i| ≥ |dot|`, for **nonzero** exact dot `cond_sum ≥ 1` and `cond_oro ≥ 2`; a finite
  `target_cond < 2` is **invalid** under the ORO convention. If the exact `dot == 0` and `S > 0`, record
  `cond_oro = cond_sum = +inf` — a valid adversarial case, but **excluded** from finite log-scale target
  matching. Each generated case records `{target_cond, measured_cond_mpfr (= cond_oro from the oracle),
  seed, scale, permutation, tier}`. The "true dot" is the **§10.2 oracle** (high precision, exponent-span
  aware, precision-doubling self-checked) — not a bare W-bit mpfr evaluation.

  *Accept:* `measured_cond_mpfr` matches `target_cond` in **log-scale** tolerance
  (`|log10(measured) − log10(target)| ≤ δ`, **default `δ = 0.25` decades** unless a test case overrides
  it), **or** the generator reports `unachievable` for that tier / length / exponent-range / target. The
  error scalings match theory — **naive** dot follows the first-order scale `~γ_n·cond_sum`, equivalently
  `~½·γ_n·cond_oro`; **`Rdot`/Dot2** follows the ORO Dot2 scaling `~u + ½·γ_n²·cond_oro` (NOT a crude
  `cond ~ 1/u` threshold) — and `vRdot` **always** encloses across all seeds, conditions, orderings, and
  tiers. (Apply the same correction to the §11.4 wording.)

- **M9 — CI / sanitizers / fast-math negative compile + runtime rounding check.** CI builds the matrix
  {native-only, MPFR-on} × {gcc, clang} × {Release, Debug}, running `ctest` incl. the separate W=53 /
  W=512 mpfr suites; an ASan+UBSan Debug build required clean; the **fast-math negative compile test**
  (`try_compile` with `-ffast-math` must **fail** via the `__FAST_MATH__` `#error`). Per §9.3, do
  **not** require an arbitrary fast-math build to fail numerically — fast-math is a contract violation
  to be **rejected**, not a numerical experiment.

  **Separately and additionally**, run a **runtime contract sanity test in ordinary (accepted-flags)
  builds**, the catch-all for dangerous flags with **no predefined macro** (`-ffp-contract=fast`,
  `-fassociative-math`, `-funsafe-math-optimizations` cannot be compile-detected). It must include:
  - a **directed-rounding test with constant-folding prevented** (use `volatile` inputs and
    `#pragma STDC FENV_ACCESS ON` where supported). Typical `double` test, with runtime (non-folded)
    values: under `FE_UPWARD`, `1.0 + 2^-53 > 1.0`; under `FE_DOWNWARD` (and nearest),
    `1.0 + 2^-53 == 1.0`. Fail if the platform/compiler ignores the rounding mode;
  - **EFT invariants** for `Rtwosum` / `Rfasttwosum` / `Rtwoproduct` (these catch FMA-contraction
    breakage that a pure rounding test may miss); the EFT self-test must cover representative **normal,
    subnormal, cancellation, large/small-scale, and signed-zero** cases;
  - an **x86 FTZ/DAZ check** where available (subnormals not flushed), matching the ROADMAP architecture
    note.

  *Accept:* CI green across the matrix; ASan/UBSan clean; `-ffast-math` compile rejected; the runtime
  contract sanity test (rounding + EFT invariants + FTZ/DAZ) green in every configuration.

- **M10 — Packaging / install / export targets.** Add `install(TARGETS ... EXPORT vmplapackTargets)`,
  header installation, `vmplapackConfig.cmake` + `vmplapackConfigVersion.cmake` so
  `find_package(vmplapack)` works. Requirements:
  - **Dependency policy:** if built with MPFR support, the installed config must
    `find_dependency(MPFR)` and `find_dependency(GMP)` (or import equivalent targets); if built
    native-only, `find_package(vmplapack)` must **not** require MPFR/GMP. `vmplapackConfig.cmake` must
    `include(CMakeFindDependencyMacro)` before any `find_dependency` call.
  - **Find-module policy:** CMake does **not** ship `FindMPFR.cmake`/`FindGMP.cmake`. The package must
    either install compatible `FindMPFR`/`FindGMP` modules **or** depend on config packages that define
    imported targets such as `MPFR::MPFR` and `GMP::GMP`. Do **not** assume CMake provides these.
  - **Usage requirements:** export the FP contract on the target — **PUBLIC** for the compiled library,
    **INTERFACE** for any interface/header-only target. Attach compiler-specific FP flags via **CMake
    generator expressions** keyed on compiler ID (`$<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:...>` etc.);
    unsupported compilers must receive equivalent strict-FP flags or be rejected by configure-time
    checks (do not put bare GCC/Clang flags in `INTERFACE_COMPILE_OPTIONS`).
  - **Include dirs:** use `$<BUILD_INTERFACE:...>` / `$<INSTALL_INTERFACE:include>`.
  - **Version policy:** use `ExactVersion` for all `0.x` releases (the API/ABI is unstable, so `0.1` and
    `0.9` must NOT be treated as compatible); switch to `SameMajorVersion` only from `1.x`, after the
    public API/ABI is declared stable. Provide a custom `vmplapackConfigVersion.cmake` if a finer rule
    is needed.
  - **Header guard is the real protection:** because a consumer can override `INTERFACE_COMPILE_OPTIONS`,
    the umbrella `#ifdef __FAST_MATH__ #error ... #endif` must be present and is what actually stops a
    fast-math consumer.
  Add a separate `tests/consumer` project that `find_package`s the installed package, links
  `vmplapack::vmplapack`, and builds a verified-dot example; it must **verify the FP flags are
  inherited** (inspect `compile_commands.json`, or a small preprocessor/compile test that fails without
  the strict-FP usage requirements). *Accept:* `cmake --install` then the consumer `find_package` + link
  + compile succeeds; the consumer inherits the FP flags (verified, not assumed); a consumer adding
  `-ffast-math` is rejected; native-only and MPFR-on packages both resolve dependencies correctly.

- **M11 — Examples / benchmark suite (frozen schema).** Add a naive-vs-compensated-vs-verified driver
  and a worked `vRresidual` example. Add a benchmark suite emitting a **frozen CSV/JSON schema** (fixed
  now for reuse in M12+) with, per case, at least:

  ```text
  schema_version,
  routine, tier, precision_bits, n, generator, seed,
  target_cond, realized_cond_oro, realized_cond_sum,
  status, enclosed,
  elapsed_ns_total, time_ns_per_item, work_items, repetitions, statistic,
  mid_error_abs, mid_error_rel,
  radius, relative_radius,
  compiler_id, compiler_version, build_type, fp_contract_flags,
  cpu, os, git_sha, rounding_backend
  ```

  Disambiguate timing (a frozen schema must): `elapsed_ns_total` = total wall time for `repetitions`,
  `time_ns_per_item` = `elapsed_ns_total / (repetitions * work_items)`, `statistic` = how repetitions
  were reduced (e.g. `median`/`mean`). `status` is required: `enclosed` (boolean) alone loses the
  `ok`/`unbounded`/`non_finite`/`invalid_input` distinction needed to compare M5/M7/M11 results. Missing-
  value conventions for routines without an enclosure (naive/compensated): `status = ok`,
  `enclosed = null`, `radius = null`, `relative_radius = null`; for `unbounded`: `radius = +inf`,
  `enclosed = false`. Represent `null` explicitly (literal `null`, or a JSON sidecar — not a blank CSV
  cell). Use `radius` (matching `Rmidrad.rad`), not `radius_width`. For a verified library,
  `radius`/`relative_radius` and "up to what condition the certificate holds" matter as much as
  throughput. Use the deterministic (M3) and random/adversarial (M8) generators. **Benchmarks are
  buildable and runnable in CI *smoke* mode only; full sweeps are opt-in and labeled separately** so CI
  stays light. *Accept:* examples build/run under all configs; the benchmark emits the frozen schema
  (incl. `schema_version`, `status`, and the disambiguated timing fields, with the null/+inf
  conventions) with per-tier timing, radius, and success-rate-vs-condition data; results documented in
  the README.

---

*After M11: RELEASE LINE A — a packaged, CI'd, benchmarked verified/accurate `sum`/`dot` library
("vMPLAPACK-core") over `float`, `double`, and MPFR, with no MPLAPACK dependency. The
verified-linear-algebra phase (M12+, where the no-MPLAPACK constraint is lifted and the verification
contract in `ROADMAP.md` applies) is specified separately when reached.*

- **M12a — Verified matvec/matmul reference.** Add the linear-algebra-layer status type
  `VerificationStatus { Verified, Unverified, InvalidInput, Unsupported }`, separate from scalar
  `Rstatus`. Add row-major point-matrix reference routines:

  ```cpp
  template <class REAL>
  VerificationStatus vRgemv_point(std::ptrdiff_t m, std::ptrdiff_t n,
                                  const REAL* A, std::ptrdiff_t lda,
                                  const REAL* x, std::ptrdiff_t incx,
                                  Rmidrad<REAL>* out);

  template <class REAL>
  VerificationStatus vRgemm_point(std::ptrdiff_t m, std::ptrdiff_t n, std::ptrdiff_t k,
                                  const REAL* A, std::ptrdiff_t lda,
                                  const REAL* B, std::ptrdiff_t ldb,
                                  Rmidrad<REAL>* C, std::ptrdiff_t ldc);
  ```

  Storage is row-major. `vRgemv_point` encloses each component `sum_j A[i,j]*x[j]`; `vRgemm_point`
  encloses each component `sum_t A[i,t]*B[t,j]`. M12a is the slow correctness baseline: own
  certification path via directed scalar component enclosures, no BLAS/MPLAPACK in the proof path, no
  M12b nearest-rounding fast enclosure. Boundary rules: negative dimensions or invalid pointers,
  strides, or leading dimensions return `InvalidInput`; zero output dimensions write nothing and return
  `Verified`; zero inner dimension writes exact zero boxes; finite regular components with
  `Rstatus::ok` yield top-level `Verified`; any `unbounded` or `non_finite` component yields top-level
  `Unverified`. `Unsupported` is reserved for future LA modes not implemented by M12a.

  *Accept:* MPFR-on tests verify `float`, `double`, and MPFR W=53/W=512 component inclusion against
  the high-precision dot oracle for cancellation-heavy matvec and matmul cases; boundary tests cover
  zero dimensions, invalid inputs, non-finite inputs, overflow-to-`unbounded`, and strided vectors.
  Native-only builds compile the float/double instantiations and keep all native tests green.

- **M12b — Verified matvec/matmul fast path.** Keep the M12a public signatures and boundary rules.
  For each component, `vRgemv_point` and `vRgemm_point` first try a nearest-rounding enclosure based
  on `vRdot_apriori`: midpoint from compensated `Rdot`, radius from the Dot2 a-priori bound using an
  upward absolute-product sum. If that a-priori path returns `unbounded`, fall back to the M12a
  directed scalar `vRdot` enclosure for that component and keep whichever result has the better scalar
  status. The certificate must not depend on an external BLAS/MPLAPACK kernel honoring directed
  rounding. `vRgemm_midrad` remains deferred.

  *Accept:* MPFR-on tests verify that a cancellation-heavy matvec/matmul component is strictly tighter
  than the directed reference while still covering the high-precision oracle interval for `float`,
  `double`, and MPFR W=53/W=512. Native tests verify that an overflowing a-priori absolute-product
  bound falls back to the directed reference instead of reporting `unbounded` for an exactly cancelling
  finite case.

- **M13 — Verified solve + nonsingularity contract.** Add `vRgesv` overloads for square point
  systems with vector and matrix right-hand sides:

  ```cpp
  template <class REAL>
  VerificationStatus vRgesv(std::ptrdiff_t n,
                            const REAL* A, std::ptrdiff_t lda,
                            const REAL* b, std::ptrdiff_t incb,
                            Rmidrad<REAL>* x, std::ptrdiff_t incx);

  template <class REAL>
  VerificationStatus vRgesv(std::ptrdiff_t n, std::ptrdiff_t nrhs,
                            const REAL* A, std::ptrdiff_t lda,
                            const REAL* B, std::ptrdiff_t ldb,
                            Rmidrad<REAL>* X, std::ptrdiff_t ldx);
  ```

  Storage is row-major. `A` is an `n x n` point matrix indexed `A[i*lda+j]`, with `lda >= n`. For
  the vector overload, `b[i*incb]` is the right-hand side and `x[i*incx]` is the output enclosure for
  solution component `i`. For the matrix overload, `B[i*ldb+r]` is right-hand side column `r` and
  `X[i*ldx+r]` is the corresponding output enclosure, with `ldb >= nrhs` and `ldx >= nrhs`. Leading
  dimensions and increments are measured in scalar elements, not bytes. Input and output ranges must
  not overlap.

  The implementation computes untrusted floating-point approximations `x~`/`X~` and an untrusted
  approximate inverse or preconditioner `R`, then certifies them. The M13 certificate is normwise in
  the infinity norm: compute an upward enclosure of `alpha = ||I - R*A||_inf`; for each RHS compute
  the residual `b - A*x~` accurately, then an upward enclosure of
  `beta = ||R*(b - A*x~)||_inf`. If `alpha < 1`, the system is certified nonsingular and every
  component of that RHS is enclosed with midpoint `x~_i` and radius `beta/(1-alpha)`. For matrix RHS,
  `alpha` is shared and `beta`/radius are per RHS column.

  Boundary rules: negative dimensions, null required pointers, `lda < n`, `incb < 1`, `incx < 1`,
  `ldb < nrhs`, or `ldx < nrhs` return `InvalidInput`. `n == 0` and, for the matrix overload,
  `nrhs == 0`, write nothing and return `Verified`. Non-finite input returns `Unverified` and writes
  `Rstatus::non_finite` boxes for all requested solution components. If `alpha >= 1`, a required
  bound overflows, or certification otherwise fails, return `Unverified` and write `unbounded` boxes;
  this is not a claim that `A` is singular. `Verified` means all requested solution components have
  `Rstatus::ok` and the true solution is included.

---

## Reference

```
T. Ogita, S. M. Rump, S. Oishi, "Accurate Sum and Dot Product",
SIAM J. Sci. Comput. 26(6), 1955–1988 (2005). DOI: 10.1137/030601818
  - Prop. 5.5: Dot2 a-priori bound; underflow term 5·n·η.
  - Cor. 5.7: relative error  u + ½·γ_n²·cond,  cond = 2·|x|ᵀ|y| / |xᵀy|  (= cond_oro).
  - Alg. 5.8 (Dot2Err): the SEPARATE pure-FP error-bound routine, gate 2·n·u < 1 — not M7.
```
