<!--
Copyright (c) 2026, NAKATA Maho
SPDX-License-Identifier: BSD-2-Clause
-->

# ROADMAP.md — vMPLAPACK

`docs/SPEC.md` specifies the prototype core, **M0–M6**. This file defines the full milestone ladder,
two release lines, an honest feasibility ranking, and the **verification contract** that governs the
whole project. A milestone is committed only when written into `docs/SPEC.md` (M7–M11 specs are included
in `docs/SPEC.md` §14 addendum; M12+ are specified when reached).

---

## Verification contract (applies to ALL milestones)

These are non-negotiable design invariants; they belong in `docs/SPEC.md` and are summarized here.

**1. Status is not a boolean.** Verified routines (the linear-algebra layer, M12+) return an explicit
status:

```cpp
enum class VerificationStatus { Verified, Unverified, InvalidInput, Unsupported };
```

`Unverified` means the routine did **not prove** the requested property with the current arithmetic,
approximation, and bounds — it does **not** mean the mathematical statement is false. (The scalar
enclosure layer M0–M11 keeps its own `Rstatus { ok, unbounded, non_finite, invalid_input }`; the two
status types serve different layers and must not be conflated.)

**2. Strict IEEE semantics required.** Native verified routines require IEEE-754 binary FP, no
fast-math, no unsafe reassociation, a controlled FMA-contraction policy, working `fenv`/`fesetround`,
subnormals not flushed (FTZ/DAZ off), and a **runtime rounding-mode sanity check**. Note: `__FAST_MATH__`
rejection catches `-ffast-math`/`-Ofast` only; `-ffp-contract=fast`, `-fassociative-math`,
`-funsafe-math-optimizations` have **no reliable predefined macro**, so the **runtime rounding-mode +
EFT self-test is the catch-all** — compile-time rejection alone is insufficient.

**3. BLAS/MPLAPACK trust boundary.** MPLAPACK/BLAS results are **approximate data only**. Certification
code must **not** rely on their internal rounding behavior (vectorization, FMA, blocking, threading, and
reassociation can all violate a requested rounding mode) unless a backend is explicitly validated as
rounding-mode safe. Correctness must come from **independently computed** residuals, rigorous norm
upper bounds, and interval/midrad enclosures.

**4. Rigorous norm primitives.** From M13 on, every norm in a bound must be a rigorous **upward** upper
bound, exposed as named primitives:

```text
vRsum_abs_upper   vRnorm_inf_upper   vRnorm_one_upper   vRmat_abs_upper   vRresidual_upper
```

---

## Milestone ladder

```
PROTOTYPE CORE (standalone, no MPLAPACK)
  M0  Skeleton + native arith
  M1  EFT, native
  M2  MPFR tier + API probe
  M3  Oracle + deterministic generators
  M4  Accurate (Rsum, Rdot)
  M5  Verified + residual (vRsum, vRdot, vRresidual)
  M6  Accurate + verified (Dot2 midpoint)

CONSOLIDATION / RELEASE PREP (still standalone, no MPLAPACK)
  M7   ORO a-priori accuracy bound
  M8   Random + adversarial high-condition generators
  M9   CI / sanitizers / fast-math negative compile + runtime rounding check
  M10  Packaging / install / export targets
  M11  Examples / benchmark suite (fixed CSV schema)
  ===> RELEASE LINE A: "vMPLAPACK-core" — standalone verified/accurate sum/dot kernels over
       float, double, MPFR. NO MPLAPACK dependency. First taggable release.

VERIFIED LINEAR ALGEBRA  (no-MPLAPACK constraint LIFTED; MPLAPACK = approximate engine only)
  M12  Verified matvec / matmul        [M12a slow directed reference; M12b fast nearest-rounding enclosure]
  M13  Verified linear solve (vRgesv)
  M14  Verified inverse / cond / det   [M14a inv; M14b cond; M14c det — det is a SEPARATE algorithm]
  ===> RELEASE LINE B: verified solve + inverse + condition bounds (the verifylss analogue). v1.0 of
       the LA goal. (Determinant M14c may trail.)
  M15  Verified positive definiteness (vRpodef)                         [stretch]
  M16  Verified Hermitian eigenvalues / invariant subspaces            [stretch]

FRONTIER (aim, NOT milestone-committed)
  vRgeev (nonsymmetric eigen), vRgesvd (SVD), vRgels (least squares), certified factor enclosures
```

---

## M7–M11 (consolidation) — summary

Full specs are in the `docs/SPEC.md` §14 addendum. These turn the M0–M6 prototype into a shippable kernel library;
all remain standalone (no MPLAPACK).

| # | Goal | Notes |
|---|------|-------|
| **M7**  | ORO a-priori accuracy bound | `vRdot_apriori` (mid=Dot2, rad=ORO bound). Rigorous form `rad ≥ outward((u·\|res\| + γ_n²·S_up + underflow_term)/(1−u))` — the `/(1−u)` accounts for `\|xᵀy\|→\|res\|`; `S_up≥Σ\|xᵢyᵢ\|` via a named upward primitive (**never** accurate `Rdot`); `underflow_term` = ORO Prop 5.5 term `5·n·η` (or a checked no-underflow precondition); for MPFR, underflow must be **checked/excluded explicitly** — do not assume the wide exponent range proves its absence. `n·u<1` gate (integer `n<2^p`); NOT Alg. 5.8 `Dot2Err` (`2n·u<1`). On gate fail return `{Rdot, +inf, unbounded}`. |
| **M8**  | Random + adversarial generators | Seeded randomized `gendot` **and** deterministic adversarial cases (cancellation, alternating signs, scale mixing, sorted/reversed/shuffled order). Record `{target_cond, measured_cond_mpfr, seed, scale, permutation, tier}`; target ≠ measured. |
| **M9**  | CI / sanitizers / contract tests | {native, MPFR}×{gcc, clang}×{Release, Debug}; ASan+UBSan clean; `-ffast-math` **try_compile must fail** (`__FAST_MATH__`); **plus a runtime rounding-mode sanity check** (FE_UPWARD changes a known expression) — the catch-all for flags that can't be macro-detected. |
| **M10** | Packaging / install / export | `install(EXPORT)`, `vmplapackConfig.cmake` + version, `find_dependency(MPFR)`/`find_dependency(GMP)` when MPFR is on, `vmplapack::vmplapack` carrying §9 flags PUBLIC/INTERFACE. The header `__FAST_MATH__` guard is the real protection (consumers can override INTERFACE flags). A `find_package` consumer smoke project. |
| **M11** | Examples / benchmark suite | Naive/compensated/verified driver; worked residual example. Benchmark records a frozen schema incl. `schema_version, status, radius, relative_radius, enclosed, realized_cond_oro/sum, time, throughput, tier, compiler, flags, cpu, os, git_sha, rounding_backend`. Use `radius` (matching `Rmidrad.rad`), not `radius_width`. **Fix the CSV schema now** for reuse in M12+. |

---

## M12–M16 (verified linear algebra) — detail

**Strategy: certification *over* MPLAPACK (the verifylss model).** Compute an approximate result with
MPLAPACK's kernels, then **certify it** with a rigorous residual bound and Rump-type fixed-point
operators. The hard targets become *result-certification*, not *verified-factorization* — you never
enclose `L,U,Q,R` themselves. Matrix representation: **midpoint–radius (midrad) matrices**, consistent
with scalar `Rmidrad<REAL>`. Tiers stay co-equal via `Rarith<REAL>`.

| # | Goal | Primitive(s) | Method (sketch) | Feasibility |
|---|------|--------------|-----------------|-------------|
| **M12a** | Verified matvec/matmul — reference | `vRgemv_point`, `vRgemm_point` | own loops, directed lower/upper passes; **no BLAS in the certification path**; correctness baseline | **Solid** |
| **M12b** | Verified matmul — fast | `vRgemm_point` (fast path), `vRgemm_midrad` (deferred) | Ozaki/Rump-style **nearest-rounding** enclosure: optimized `Rgemm` for the approximate midpoint, error bound from a-priori analysis — **must not assume BLAS honors directed rounding** | **Solid, with care** |
| **M13** | Verified solve + nonsingularity | `vRgesv` (Ax=b, AX=B) | approx x̃ and approx inverse/preconditioner R (FP); enclose α=‖I−RA‖ (upward); if α<1 then ‖x−x̃‖ ≤ ‖R(b−Ax̃)‖/(1−α); **residual `b−Ax̃` computed accurately** (M12 / compensated dot — essential for ill-conditioned A). α≥1 ⇒ `Unverified`, **not** "A singular". Consider componentwise bounds (normwise can be pessimistic). | **Solid** |
| **M14a** | Verified inverse | `vRgeinv` | solve `AX=I`; return a midrad enclosure of A⁻¹ | **Solid** |
| **M14b** | Certified condition bounds | `vRgecon` | return certified bounds with **explicitly specified meaning** (e.g. upper bound on ‖A‖·‖A⁻¹‖ and/or lower bound on rcond); do not claim an exact condition enclosure | **Solid-ish** |
| **M14c** | Verified determinant | `vRgedet` | **separate algorithm**; current reference certificate directly encloses the Leibniz expansion for `n <= 8`, returning `unbounded` above that. A scalable Rump-style bound can replace the internals later. **Not** a byproduct of M13. | **Reference baseline** |
| **M15** | Verified positive definiteness | `vRpodef` (SPD certificate + λ_min lower bound) | FP Cholesky-based **candidate** for A−cI, then verify **Rump's sufficient criterion including all rounding errors** (Rump 2006). *Completion of FP Cholesky alone is NOT a certificate.* If certified ⇒ A is SPD and λ_min(A) > c. Fix the input-symmetry policy (assume exact symmetric storage, or check + symmetrize with an explicit enclosure). | **Stretch** |
| **M16** | Verified Hermitian eigenvalues / invariant subspaces | `vRsyev`/`vRheev`, `vReigpair` | (1) eigenvalue intervals via a **verified Hermitian inertia routine** (signs of a verified LDLᵀ/LDL* of A−σI with rigorous error control) + bisection — A−σI is generally **indefinite**, so Cholesky does **not** give full inertia; the M15 PD certificate only serves threshold tests (all eigenvalues above/below σ), not counting; (2) simple eigenpair residual bounds (Rayleigh quotient, Weyl/Temple); (3) **clustered invariant-subspace** enclosures via Hermitian-specific methods (Rump–Lange 2023). Bauer–Fike belongs to the nonsymmetric frontier, not here. | **Stretch** |

**Release line B = M12 + M13 + M14a/b** (verified solve, inverse, condition). v1.0 of the LA goal;
`vRgedet` (M14c) is separate and currently reference-quality.

---

## Frontier — aim, but NOT milestone-committed

Litmus test before promoting: *can it be expressed as certifying an MPLAPACK-computed result via a
residual + a fixed-point operator with a computable contraction factor?* If yes → likely tractable; if
it needs enclosing internal factors → defer.

- **Nonsymmetric eigenvalues** (`vRgeev`): clustered/defective/non-normal are hard. Simple eigenvalues
  reachable (Rump 2001); all eigenvalues/eigenvectors incl. clusters: Rump 2022.
- **SVD** (`vRgesvd`): reduction to a symmetric eigenproblem of AᵀA is mathematically valid but **squares
  the condition number** and should **not** be the primary verified algorithm. Prefer block-Hermitian
  formulations or residual-based singular-pair certification (Rump–Lange 2023).
- **Certified factor enclosures** (`L,U`/`Q,R` themselves): rarely needed, harder than result
  certification — out of scope unless a use case demands it.
- **Ill-conditioned least squares** (`vRgels`): augmented-system certification possible but delicate.

---

## Architectural notes

- **No-MPLAPACK boundary:** M0–M11 are standalone. **From M12 the constraint is lifted** — MPLAPACK
  approximate kernels (`Rgemm`, `Rgetrf`, `Rgetri`, `Rpotrf`, `Rsyev`/`Rheev`, …) are the inner engine;
  vMPLAPACK is the certifying wrapper. Their outputs are **untrusted approximate data** (contract §3).
- **Midrad matrix type** parallel to `Rmidrad<REAL>`; the **verified matmul** is the cornerstone. Start
  with point × point enclosures (`vRgemm_point`); midrad × midrad is deferred.
- The M13 residual `b − Ax̃` must be computed **accurately** (compensated / higher precision) — for
  ill-conditioned systems a naively-rounded residual is meaningless. This is the structural reason the
  accurate/verified dot kernels (M0–M11) precede the LA layer.
- **x86 specifics:** FTZ/DAZ, FMA contraction, and x87-vs-SSE differences all affect soundness; the
  rounding contract (§9) and the runtime sanity check must cover them.

---

## References (DOIs verified)

```
Core (M0–M11):
- T. Ogita, S. M. Rump, S. Oishi, "Accurate Sum and Dot Product",
  SIAM J. Sci. Comput. 26(6), 1955–1988 (2005).            DOI: 10.1137/030601818
- N. J. Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., SIAM (2002).
  ISBN: 0-89871-521-0 / 978-0-89871-521-7.                 DOI: 10.1137/1.9780898718027

Verified matmul (M12):
- K. Ozaki, T. Ogita, S. M. Rump, S. Oishi, "Fast algorithms for floating-point interval
  matrix multiplication", J. Comput. Appl. Math. 236(7), 1795–1814 (2012).
                                                            DOI: 10.1016/j.cam.2011.10.011
- S. M. Rump, "Fast interval matrix multiplication",
  Numerical Algorithms 61, 1–34 (2012).                    DOI: 10.1007/s11075-011-9524-z
- S. M. Rump, "Fast and Parallel Interval Arithmetic",
  BIT Numerical Mathematics 39, 534–554 (1999).            DOI: 10.1023/A:1022374804152

Verified solve / inverse / determinant (M13–M14):
- S. Oishi, S. M. Rump, "Fast verification of solutions of matrix equations",
  Numerische Mathematik 90(4), 755–773 (2002).             DOI: 10.1007/s002110100310
- S. M. Rump, "Verified bounds for the determinant of real or complex point or interval
  matrices", J. Comput. Appl. Math. 372, 112610 (2020).    DOI: 10.1016/j.cam.2019.112610

Verified definiteness / eigen (M15–M16, frontier):
- S. M. Rump, "Verification of Positive Definiteness",
  BIT Numerical Mathematics 46, 433–452 (2006).            DOI: 10.1007/s10543-006-0056-1
- S. M. Rump, "Computational Error Bounds for Multiple or Nearly Multiple Eigenvalues",
  Linear Algebra Appl. 324, 209–226 (2001).                DOI: 10.1016/S0024-3795(00)00279-2
- S. M. Rump, "Verified Error Bounds for All Eigenvalues and Eigenvectors of a Matrix",
  SIAM J. Matrix Anal. Appl. 43(4), 1736–1754 (2022).      DOI: 10.1137/21M1451440
- S. M. Rump, M. Lange, "Fast computation of error bounds for all eigenpairs of a Hermitian
  and all singular pairs of a rectangular matrix with emphasis on eigen- and singular value
  clusters", J. Comput. Appl. Math. 434, 115332 (2023).   DOI: 10.1016/j.cam.2023.115332

Survey (all phases):
- S. M. Rump, "Verification methods: Rigorous results using floating-point arithmetic",
  Acta Numerica 19, 287–449 (2010).                        DOI: 10.1017/S096249291000005X
```
*(DOIs cross-checked against publisher/author records, 2026-06-23; re-confirm before public release.)*
