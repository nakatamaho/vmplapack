<!--
Copyright (c) 2026, NAKATA Maho
SPDX-License-Identifier: BSD-2-Clause
-->

# SKILLS.md — vMPLAPACK implementation playbook

Concrete, reusable patterns for the parts of vMPLAPACK that are easy to get *plausibly* right and
*actually* wrong. Read before implementing kernels. Each skill gives the **rule**, the **correct
pattern**, and the **anti-pattern**. The authoritative spec is `docs/SPEC.md`; this file is the
"how", not the "what".

Overarching principle: **a wrong result that looks bounded is worse than an honestly unbounded
one.** Prefer `unbounded` to a false certificate; prefer stopping to guessing.

---

## Skill 1 — Error-free transformations (EFT)

**Rule.** EFTs must be **alias-safe** (an output may be the same variable as an input) and are
**exact only under preconditions**.

**Correct pattern — pass inputs by value:**

```cpp
// a + b = s + e exactly (Knuth TwoSum), s = fl(a+b). Branch-free, alias-safe.
template <class REAL>
void Rtwosum(REAL a, REAL b, REAL& s, REAL& e) {
    s = a + b;
    REAL z  = s - a;
    REAL t1 = s - z;
    REAL t2 = a - t1;
    REAL t3 = b - z;
    e = t2 + t3;
}

// a * b = p + e exactly (under preconditions), p = fl(a*b). Alias-safe.
template <class REAL>
void Rtwoproduct(REAL a, REAL b, REAL& p, REAL& e) {
    p = a * b;
    REAL neg_p = -p;
    e = Rarith<REAL>::fma(a, b, neg_p);   // correctly-rounded fma: std::fma / fmaf / mpfrxx::fma
}
```

The reductions call these aliased: `Rtwosum(s, x[i], s, e)`, `Rtwosum(p, h, p, q)`. **This must
work.** By-value inputs make it work; the MPFR copy cost is accepted.

**Exactness preconditions (do not over-claim):**
- `Rtwosum`: finite inputs, no overflow in `a+b`, subnormals not flushed, round-to-nearest.
- `Rfasttwosum`: additionally `abs(a) >= abs(b)`.
- `Rtwoproduct`: additionally the residual `a*b - fl(a*b)` must be representable — i.e. the product
  must **not overflow AND not underflow**. Under product underflow, `e` is the *correctly-rounded*
  residual but **not exact**.

**Anti-patterns.**
- ❌ `void Rtwosum(const REAL& a, ..., REAL& s, ...)` then read `a` after writing `s` — when `s`
  aliases `a`, the first `s = a + b` corrupts `a`; every later `s - a` is wrong.
- ❌ Asserting `a*b == p+e` for underflowing or overflowing products. Test those as **robustness**
  cases (correct status / valid enclosure), not as exact-TwoProduct cases.
- ❌ Dekker `Rsplit` with `factor*a` near the format max — it overflows. Split factor is
  `2^ceil(p/2)+1`: float `0x1p12f+1`, double `0x1p27+1`, mpfr `2^ceil(W/2)+1`.

---

## Skill 2 — The rounding / floating-point contract

**Rule.** Verification is valid only if the compiler/runtime does not rewrite the arithmetic and
the rounding mode is actually in force. This is **load-bearing**.

**Correct pattern — per-tier rounding via `Rarith<REAL>`, RAII, non-copyable:**

```cpp
// native (float/double): fesetround; copy/move deleted so we never restore twice
struct round_up {
    int saved_;
    round_up() : saved_(std::fegetround()) { std::fesetround(FE_UPWARD); }
    ~round_up() { std::fesetround(saved_); }
    round_up(const round_up&) = delete;  round_up& operator=(const round_up&) = delete;
    round_up(round_up&&) = delete;       round_up& operator=(round_up&&) = delete;
};
// round_down is identical with FE_DOWNWARD.

// mpfr: wrap the library RAII scope (already non-copyable)
struct round_up   { mpfrxx::rounding_mode_scope s_{MPFR_RNDU}; };
struct round_down { mpfrxx::rounding_mode_scope s_{MPFR_RNDD}; };
```

**Build/contract checklist (native TUs):**
- Flags: `-ffp-contract=off -fno-fast-math -frounding-math` (MSVC `/fp:strict`), exposed
  PUBLIC/INTERFACE so consumers inherit them.
- In each such TU: `#pragma STDC FENV_ACCESS ON` and `#pragma STDC FP_CONTRACT OFF`.
- Umbrella header: `#if defined(__FAST_MATH__) #error ... #endif` — **compile-time rejection is the
  primary safeguard.**
- `static_assert(FLT_EVAL_METHOD == 0, ...)`; on 32-bit x86 add `-mfpmath=sse -msse2`.
- **FTZ/DAZ off** (subnormals alive).

**Runtime self-test (mandatory, per tier):** confirm (a) mode restored to nearest after a scope and
on exception; (b) crafted EFT reconstruction exact; (c) **subnormals survive** AND entering
`round_up`/`round_down` **actually changes** an ordinary `a + b` and `a * b` at runtime; (d) directed
passes bracket a known hard case.

**Anti-patterns.**
- ❌ Trusting `-frounding-math` alone (GCC documents it as experimental). Use flag **+** pragma **+**
  runtime test.
- ❌ Calling EFT/accurate kernels inside a directed scope — they require round-to-nearest.
- ❌ Letting `acc + x[i]*y[i]` contract into an FMA in a verified pass (breaks the per-op rounding
  model). Keep `-ffp-contract=off`; use plain product + add in verified passes (FMA only in
  `Rtwoproduct`).
- ❌ Requiring an arbitrary `-ffast-math` build to fail *numerically* — it may not; rely on the
  compile-time `__FAST_MATH__` rejection.

---

## Skill 3 — Verified enclosure construction

**Rule.** Two directed passes give `true ∈ [inf, sup]`. Convert to a midrad enclosure with
**overflow-safe** arithmetic and a **single source of truth** for NaN/Inf classification.

**Correct pattern — directed two-pass + helpers:**

```cpp
// upper bound: every product and add rounded toward +inf
{ Rarith<REAL>::round_up   g; acc = zero(); for i: { prod = x[ix]*y[iy]; acc = acc + prod; } sup = acc; }
// lower bound: every product and add rounded toward -inf
{ Rarith<REAL>::round_down g; acc = zero(); for i: { prod = x[ix]*y[iy]; acc = acc + prod; } inf = acc; }

REAL mid = Rmidpoint(inf, sup);            // lo/2 + hi/2 — never overflows
return Rmake_midrad(inf, sup, mid);        // classifies bounds -> ok / unbounded
```

```cpp
template <class REAL> REAL Rmidpoint(REAL lo, REAL hi) {
    REAL h = Rarith<REAL>::half(); REAL a = lo*h; REAL b = hi*h; REAL m = a + b; return m;
}
template <class REAL> REAL Rupward_abs_diff(REAL a, REAL b) {        // |a-b| rounded UP
    typename Rarith<REAL>::round_up g; REAL d = (a>=b)?(a-b):(b-a); return d;
}
template <class REAL> Rmidrad<REAL> Rmake_midrad(REAL lo, REAL hi, REAL mid) {
    using A = Rarith<REAL>;
    if (!A::is_finite(lo) || !A::is_finite(hi)) return { A::zero(), A::infinity(), Rstatus::unbounded };
    if (hi < lo) { assert(false && "lo>hi"); return { A::zero(), A::infinity(), Rstatus::unbounded }; }
    REAL m = mid;
    if (!A::is_finite(m)) m = Rmidpoint(lo, hi);   // keep the finite enclosure
    REAL r1 = Rupward_abs_diff(m, lo), r2 = Rupward_abs_diff(hi, m);
    REAL rad = (r1 > r2) ? r1 : r2;
    if (!A::is_finite(rad)) return { A::zero(), A::infinity(), Rstatus::unbounded };
    return { m, rad, Rstatus::ok };
}
```

`Rmake_midrad` is rigorous even when `mid` is outside `[lo, hi]` — one upward distance covers the
whole interval, so `[mid-rad, mid+rad] ⊇ [lo, hi] ∋ true`.

**NaN/Inf classification — the single rule (`SPEC §6.3`):**
1. **First scan inputs.** Any NaN/Inf **input** → `non_finite`, **no certificate** (`is_finite` is
   the only predicate you need — false for both NaN and Inf).
2. After inputs are finite, the true value is finite, so **any** non-finite endpoint/mid/radius is
   overflow → `unbounded` with `mid=0`, `rad=+inf` (the valid, useless `[-inf,+inf]`).
3. `Rmake_midrad` is **not** an input classifier; it only handles already-finite-input cases.

**M6 (accurate + verified):** `mid = Rdot(x,y)`, then `Rmake_midrad(inf, sup, mid)`. The radius stays
rigorous regardless of where `mid` lands; if `mid` overflows while `[inf,sup]` is finite, the helper
falls back to `Rmidpoint` and keeps the finite enclosure.

**Anti-patterns.**
- ❌ `mid = (inf + sup) / 2` — overflows when both bounds are large. Use `lo/2 + hi/2`.
- ❌ Classifying a NaN endpoint as `non_finite` — that contradicts the helper and the spec. Finite
  inputs ⇒ overflow ⇒ `unbounded`.
- ❌ Returning `ok` without checking `lo <= hi` — a contract breach would slip through.
- ❌ Discarding a finite `[inf,sup]` just because the accurate midpoint failed.

---

## Skill 4 — MPFR tier via gmpfrxx_mkII

**Rule.** mpfr is a **first-class tier**, at a **fixed** precision `W`. Respect expression-template
evaluation, fixed-`W` semantics, and signed exponents.

**Correct patterns.**

```cpp
// Materialize every elementary op (operators return EXPRESSION NODES, not mpfr_class):
mpfrxx::mpfr_class prod = xi * yi;   // force evaluation into a value
acc = acc + prod;                    // force again — NEVER `auto r = xi*yi;`

// Correctly-rounded TwoProduct on mpfr:
e = mpfrxx::fma(a, b, neg_p);

// 2^-W with a SIGNED exponent (precision_bits() returns signed `long`):
mpfrxx::mpfr_class u = mpfrxx::mpfr_class::with_precision(W);
mpfr_exp_t ee = -static_cast<mpfr_exp_t>(W);     // signed before negation
mpfr_set_ui_2exp(u.mpfr_data(), 1, ee, MPFR_RNDN);

// Exact powers for test data (never via double):
mpfrxx::mpfr_class p2j = mpfrxx::mpfr_class::with_precision(W);
mpfr_set_ui_2exp(p2j.mpfr_data(), 1, j, MPFR_RNDN);   // = 2^j exactly
```

**Fixed-`W` discipline (`SPEC §10.1`):** MPFR's default precision only affects *subsequently*
initialized objects. So: set `mpfrxx::set_default_precision_bits(W)` **before** constructing any
input; every input must have precision `W` (assert `mpfr_get_prec(x.mpfr_data()) == W`); **never**
change `W` after any `mpfr_class` exists. Run `W=53` and `W=512` as **separate** test processes.

**API surface used** (verified against gmpfrxx_mkII source; **M2 probes it first**):
`mpfrxx::mpfr_class`, `with_precision(prec[, value])`, `default_prec()`,
`set_default_precision_bits()`, `mpfr_data()`/`get_mpfr_t()`, `mpfrxx::fma`, `mpfrxx::abs`,
`mpfrxx::rounding_mode_scope(MPFR_RNDU/RNDD)`. Do **not** enable MPFR FMA auto-fusion; call
`mpfrxx::fma` explicitly and use plain product+add in verified passes.

**Anti-patterns.**
- ❌ `auto r = a + b;` for mpfr — captures a node referencing temporaries; may dangle, and bypasses
  per-op rounding.
- ❌ `-precision_bits()` when the type is unsigned — wraps to a huge positive exponent.
- ❌ Constructing `2^j` via `std::ldexp(double(1), j)` for the mpfr tier — overflows/loses precision
  through `double`.
- ❌ Changing `W` mid-process or mixing `W=53` and `W=512` objects in one run.

---

## Skill 5 — The oracle and condition-aware testing

**Rule.** The oracle is a **naive high-precision interval**, self-certified by precision doubling.
Accurate tests must be **condition-aware**; verified tests assert **interval inclusion**.

**Correct patterns.**

```text
Oracle (tests only): naive (no EFT) dot at precision P, returned as an INTERVAL:
    ref_lo = naive dot at P, MPFR_RNDD     // <= true
    ref_hi = naive dot at P, MPFR_RNDU     // >= true
Start P = max(512, 4*p + 64); the RIGOROUS check is doubling: also compute at 2P and require the two
intervals to agree far below u_tier. Widen P until it passes. Exact input widening: with_precision(P, d)
for native; mpfr_class(x.mpfr_data(), P) for the mpfr tier (P >= W) — never move across precisions.
```

```text
Accurate test (condition-aware; NOT |got-ref| <= C*u*|ref|):
    S_hi = oracle upper bound of Σ|x_i y_i| via an UPWARD-directed high-precision pass
    ref  = midpoint of [ref_lo, ref_hi]
    require |got - ref| <= u*|ref| + C*u^2*S_hi + floor      // C ~ 4 (test heuristic)
    floor = a few * (ref_hi - ref_lo)                        // oracle width; tier-agnostic
This passes Dot2 yet fails a non-compensated dot (error ~ u*S ≫ u^2*S when S ≫ |ref|).
```

```text
Verified inclusion test:
    require status in {ok, unbounded};  rad >= 0
    require (mid - rad) <= ref_lo  AND  ref_hi <= (mid + rad)    // covers the whole oracle interval
Overflow inputs -> require unbounded; non-finite inputs -> require non_finite. Never a false ok.
```

**Anti-patterns.**
- ❌ `|got - ref| <= C*u*|ref|` — breaks on cancellation / exact-zero `ref` (right side → 0).
- ❌ Using accurate `Rdot(|x|,|y|)` for `S` in a bound — it can under-estimate. Use an upward pass.
- ❌ A native "smallest-normal" floor for the mpfr tier — there is no such thing; use the
  oracle-width floor.
- ❌ Inclusion against a single high-precision **point** — the point itself carries oracle rounding;
  cover the **interval** `[ref_lo, ref_hi]`.

---

## Skill 6 — Ill-conditioned test data (`Rgendot.h`)

**Rule.** Generated dot cases must have reproducible metadata and exactly constructed tier inputs.
Condition-number targets use the ORO convention, not an implicit relative-error threshold.

```text
cond_oro = 2 * sum_i |x_i*y_i| / |x'y|     # ORO 2005 convention
cond_sum =     sum_i |x_i*y_i| / |x'y|     # sum convention
```

For nonzero exact dots, `cond_oro >= 2`. A finite target below 2 is invalid. For a nontrivial exact
zero dot, record `cond_oro = cond_sum = +inf` and do not try to match a finite log-scale target.

**Correct patterns.**

```text
Family A (alternating):  x_i = 1; y chosen so exact dot = delta.
Family B (two-term):     x = [1, 1], y = [1, y2], y2 exactly representable; exact dot = 1 + y2.
Family C (exponent-cancellation):
    x = [2^j, s, 2^j], y = [1, 1, -1]
    Place s between the large terms so naive sequential accumulation loses s when j >= W.
M8 adversarial families:
    heavy cancellation, alternating signs, huge/small scale mixing, exact nontrivial cancellation.
Ordering variants:
    generated, sorted, reversed, shuffled variants of the same multiset.
```

The seeded high-condition generator targets a finite `cond_oro` and records:

```text
target_cond, measured_cond_mpfr, measured cond_sum, seed, scale, permutation, tier, status
```

`measured_cond_mpfr` is computed from the §10.2 oracle interval and an upward high-precision absolute
term sum, not from a bare W-bit mpfr dot. Acceptance is log-scale:

```text
|log10(measured_cond_mpfr) - log10(target_cond)| <= 0.25 decades
```

or the generator reports `unachievable` for that tier/length/exponent range/target.

**Exact-construction rule.** Native powers use in-range `std::ldexp` in the target type. MPFR powers use
`mpfr_set_ui_2exp` at the current fixed `W`. Never construct generator input values through `double`.

**Scaling checks.** Test the variables in the theory, not a crude `cond ~ 1/u` cutoff:

```text
naive dot scale:  gamma_n * cond_sum = 0.5 * gamma_n * cond_oro
Dot2/Rdot scale:  u + 0.5 * gamma_n^2 * cond_oro
```

Across all seeds, conditions, deterministic families, and orderings, `vRdot` must enclose the oracle
interval. `Rdot` may lose accuracy; verified inclusion must not fail.

**Anti-patterns.**
- ❌ Reporting only one unnamed `cond` value — the factor-of-two convention matters.
- ❌ Matching finite targets by ordinary relative error rather than log-scale tolerance.
- ❌ Computing `S` or measured condition from accurate `Rdot`; use an upward absolute-term pass.
- ❌ Letting shuffled/sorted/reversed variants change the multiset rather than only the order.

---

## Skill 7 — Debugging heuristics & final anti-pattern checklist

**Heuristics.**
- **"~2^-53 error at high precision" ⇒ a hidden `double` detour.** If an mpfr@512 (or @W) result is
  accurate only to ~`2^-53` while `u = 2^-W` is far smaller, suspect an accidental conversion to
  `double` — e.g. a missing overload causing ADL to resolve a call to `std::` via an implicit
  `operator double()`. Diff the mpfr path against the native path; look for a bare `a + b`/`a * b`
  that silently truncated to 53 bits.
- **Inclusion fails only at extreme exponent gaps ⇒ contraction or FTZ.** Re-check `-ffp-contract=off`
  and that subnormals are alive; an unexpected FMA or a flushed subnormal breaks a directed bound.
- **`unbounded` where you expected `ok` ⇒ a real overflow in the passes**, not a bug in the helper —
  confirm with the magnitudes before "fixing" the classifier.
- **EFT test passes in-tier but the kernel is still inaccurate ⇒ you tested `hi+lo` in-tier**
  (re-rounded). Re-verify exactness in MPFR.

**One-line checklist before opening a PR.**
- [ ] EFTs by value; aliased calls tested.
- [ ] Verified passes: plain product+add, directed scope per pass, materialized.
- [ ] Midpoint overflow-safe; NaN/Inf classified per `SPEC §6.3`; `lo<=hi` checked.
- [ ] FP flags + pragmas + `__FAST_MATH__` guard present; runtime contract test green.
- [ ] mpfr: fixed `W`, inputs at `W`, signed `2^-W`, exact `2^j`, materialized ops.
- [ ] Accurate test condition-aware with an **upper-bound** `S`; verified test covers the oracle
      **interval**.
- [ ] No accurate `Rdot` used as a bound; no in-tier `hi+lo` exactness test.
- [ ] Boundary cases (`n==0`, `n<0`, null, `vRresidual m==0/n==0`) per `SPEC §8`.
- [ ] No FP flag relaxed to silence a warning; no API invented; no frozen contract changed without a
      flagged report.
