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

## References

- T. Ogita, S. M. Rump, S. Oishi, "Accurate Sum and Dot Product", SIAM J. Sci. Comput. 26(6),
  1955-1988, 2005. DOI: 10.1137/030601818.
- `docs/SPEC.md`, M7: ORO a-priori accuracy bound.
