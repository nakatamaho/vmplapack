<!--
Copyright (c) 2026, NAKATA Maho
SPDX-License-Identifier: BSD-2-Clause
-->

# AGENTS.md — vMPLAPACK

Operational contract for coding agents (Codex et al.) working in this repository.

**Read this file fully, then read `SKILLS.md` before writing any numerical kernel.**
The authoritative, complete requirements live in `docs/SPEC.md`. This file is the
quick-reference for workflow, commands, and the non-negotiable rules; `SKILLS.md` is
the implementation playbook for the tricky numerical/C++ parts. When this file and
`docs/SPEC.md` disagree, `docs/SPEC.md` wins — and you must stop and flag the conflict.

---

## What this is

vMPLAPACK is a **standalone C++17 prototype** of **accurate** and **verified** `sum`/`dot`
kernels, in the MPLAPACK idiom (`template <class REAL>`). Three scalar tiers are **co-equal**,
unified by a single `Rarith<REAL>` policy:

```
REAL = float                   // IEEE-754 binary32, u = 2^-24
REAL = double                  // IEEE-754 binary64, u = 2^-53
REAL = mpfrxx::mpfr_class       // MPFR via gmpfrxx_mkII, FIXED working precision W
```

It is **NOT** an MPLAPACK patch and must not depend on the MPLAPACK source tree.

- **Accurate** (`Rsum`, `Rdot`): better approximation (compensated), returns a bare `REAL`. No certificate.
- **Verified** (`vRsum`, `vRdot`, `vRresidual`): rigorous midpoint–radius enclosure
  `true ∈ [mid-rad, mid+rad]` plus an `Rstatus`. "Verified" means **interval inclusion**, never
  "close to MPFR".

---

## Golden rules (non-negotiable)

1. **Soundness over tightness; correctness over cleverness.**
2. **Milestone order is law.** Implement strictly M0 → M12 (see `docs/SPEC.md` §14 and `ROADMAP.md`). **One PR per
   milestone.** Do not start a milestone until the previous one's acceptance tests pass.
3. **Frozen contracts.** Do not change the public signatures, `Rmidrad`/`Rstatus`, the
   `Rarith<REAL>` interface, the EFT preconditions (`SPEC §6`), the boundary rules (`SPEC §8`), or
   the rounding/FP contract (`SPEC §9`) **without stopping and reporting**.
4. **The rounding/FP contract is load-bearing for the word "verified."** Never relax it to make a
   build pass, satisfy a platform, or silence a warning.
5. **When blocked or when a requirement looks wrong: STOP and report.** Do not improvise around a
   soundness rule, do not invent an API, do not weaken a guarantee to make tests green.
6. All code, identifiers, and comments in **English**.

---

## Environment / dependencies

The MPFR tier and the test oracle need GMP + MPFR (+ MPC if you include the combined
`gmpfrxx_mkII.h`) and the **gmpfrxx_mkII** headers.

```bash
# Debian/Ubuntu
sudo apt-get update && sudo apt-get install -y cmake g++ libgmp-dev libmpfr-dev libmpc-dev
```

Provide gmpfrxx_mkII headers either as a git submodule at `extern/gmpfrxx_mkII` or via
`-DGMPFRXX_MKII_INCLUDE_DIR=/path/to/gmpfrxx_mkII/include`. **Do not vendor or modify
gmpfrxx_mkII sources.** Its public API (`mpfrxx::mpfr_class`, `mpfrxx::fma`,
`mpfrxx::rounding_mode_scope`, `mpfr_data()`, `with_precision`, `default_prec`,
`set_default_precision_bits`) is verified against its source, but **M2 must begin with a compile
probe** that confirms these symbols exist and stop with the exact missing symbol if not.

---

## Commands

```bash
# Full build (both tiers + oracle)
cmake -S . -B build -DVMPLAPACK_ENABLE_MPFR=ON -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j
ctest --test-dir build --output-on-failure        # runs float, double, mpfr@53, mpfr@512

# Run only one mpfr precision (the W=53 / W=512 suites are separate CTest cases)
MPFRXX_DEFAULT_PRECISION_BITS=53  ctest --test-dir build -R 'mpfr.*w53'  --output-on-failure
MPFRXX_DEFAULT_PRECISION_BITS=512 ctest --test-dir build -R 'mpfr.*w512' --output-on-failure

# Native-only build (no GMP/MPFR; oracle/MPFR-exactness tests are skipped by design)
cmake -S . -B build-native -DVMPLAPACK_ENABLE_MPFR=OFF
cmake --build build-native -j
ctest --test-dir build-native --output-on-failure

# Optional: force the Dekker TwoProduct path instead of FMA
cmake -S . -B build -DVMPLAPACK_USE_DEKKER_TWOPRODUCT=ON
```

`-Wall -Wextra -Werror` is expected. **Never** disable the strict FP flags
(`-ffp-contract=off -fno-fast-math -frounding-math`, `FLT_EVAL_METHOD == 0`) to clear a warning or
an error; fix the code instead.

---

## Repository layout

```
CMakeLists.txt
cmake/RoundingControl.cmake        # strict FP flags, exposed PUBLIC/INTERFACE
include/vmplapack/
  vmplapack.h                      # umbrella; carries the __FAST_MATH__ #error guard
  vmpblas.h                        # public decls + extern template decls
  vmplapack_utils.h                # Rstatus, Rstatus_rank, Rmidrad, midrad helpers
  vmplapack_arith.h                # Rarith<REAL>
  vmplapack_eft.h                  # Rtwosum / Rfasttwosum / Rtwoproduct / Rsplit
src/                               # one routine per .cpp; explicit template instantiation
  Rsum.cpp Rdot.cpp vRsum.cpp vRdot.cpp vRresidual.cpp vRgemv_point.cpp vRgemm_point.cpp
tests/                            # Rdot_oracle.h, Rgendot.h, test_*.cpp
examples/                         # example_vRdot.cpp, example_vRresidual.cpp
docs/SPEC.md                       # the authoritative specification
AGENTS.md  SKILLS.md  README.md
```

One public routine per `src/*.cpp`. Template definitions live in headers; first-class tiers get
explicit instantiation in `src/` + `extern template` in `vmpblas.h` (guard the mpfr instantiation
with `VMPLAPACK_ENABLE_MPFR`).

---

## Workflow per milestone

1. Read the milestone in `docs/SPEC.md §14` and the relevant `SKILLS.md` sections.
2. Implement only that milestone's scope.
3. Add tests; **each test comment states the failure it is designed to catch.** Verified tests
   assert **inclusion of the oracle interval `[ref_lo, ref_hi]`**, not closeness to a point.
4. Update `ROUTINES.md` for the milestone: document every product-level public routine added or changed,
   plus the important helper/subordinate routines that make the public routine work. Cover purpose,
   signatures, status behavior, runtime/oracle requirements, and how users should interpret results.
5. Run the milestone's acceptance commands; all green before proceeding.
6. Open **one PR** for the milestone.

If a milestone's acceptance cannot be met without violating a golden rule, **stop and report the
exact conflict** rather than weakening the rule.

---

## Style & conventions

- C++17. Naming: accurate `Rsum`/`Rdot`, verified `vRsum`/`vRdot`/`vRresidual`, EFT `Rtwosum`/
  `Rfasttwosum`/`Rtwoproduct`/`Rsplit`. `R` is a "real, templated-over-`REAL`" family marker — never
  put `float`/`double`/`mpfr` in a public name.
- **Materialization rule:** one rounded elementary op per statement; never bind an arithmetic
  result with `auto` (expression templates). See `SKILLS.md`.
- Keep `float`/`double`/`mpfr_class` strictly co-equal: one generic template per kernel; **all**
  tier differences live in `Rarith<REAL>`.
- Examples and demos must be tier-complete: every example must exercise `float`, `double`, and
  `mpfrxx::mpfr_class` at fixed `W=512`. Before the MPFR tier exists (M0/M1), any temporary example
  must be clearly treated as incomplete and updated to include MPFR@512 as soon as M2 lands.
- Every non-generated repository file must start with a BSD-2-Clause file header using the file's
  comment syntax: `Copyright (c) 2026, NAKATA Maho` and `SPDX-License-Identifier: BSD-2-Clause`.
  Do not add that header to `LICENSE` itself or to empty `.gitkeep` marker files.

## Commit / PR conventions

- One PR per milestone; title `Mk: <concise summary>` (e.g. `M4: accurate Rsum/Rdot (Sum2/Dot2)`).
- English commit messages; reference the milestone acceptance criteria; no unrelated changes.
- PRs that touch a frozen contract (rule 3) must say so explicitly and explain why.

---

## Do NOT

- Do **not** use `-ffast-math`, unsafe reassociation, or FP contraction; do not edit the FP flags.
- Do **not** call accurate/EFT kernels from inside a directed-rounding scope (they require
  round-to-nearest).
- Do **not** use accurate `Rdot` as an upper bound (it can under-estimate) — use an upward-directed
  pass or a `vRdot` upper endpoint.
- Do **not** test EFT exactness with same-tier `hi + lo` (it re-rounds) — compare in MPFR.
- Do **not** write `-precision_bits()` if `precision_bits()` could be unsigned — cast to a signed
  `mpfr_exp_t` first.
- Do **not** claim a finite verified certificate for an overflowed computation — return `unbounded`.
- Do **not** introduce a virtual per-scalar-op backend or per-type kernels.
- Do **not** vendor or modify gmpfrxx_mkII / MPLAPACK sources.
