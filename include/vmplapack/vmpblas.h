// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <cassert>
#include <cstddef>

#include <vmplapack/vmplapack_arith.h>
#include <vmplapack/vmplapack_eft.h>
#include <vmplapack/vmplapack_utils.h>

namespace vmplapack {

template <class REAL>
REAL Rsum(std::ptrdiff_t n, const REAL* x, std::ptrdiff_t incx = 1);

template <class REAL>
REAL Rdot(std::ptrdiff_t n,
          const REAL* x,
          std::ptrdiff_t incx,
          const REAL* y,
          std::ptrdiff_t incy);

template <class REAL>
Rmidrad<REAL> vRsum(std::ptrdiff_t n, const REAL* x, std::ptrdiff_t incx = 1);

template <class REAL>
Rmidrad<REAL> vRdot(std::ptrdiff_t n,
                    const REAL* x,
                    std::ptrdiff_t incx,
                    const REAL* y,
                    std::ptrdiff_t incy);

template <class REAL>
Rstatus vRresidual(std::ptrdiff_t m,
                   std::ptrdiff_t n,
                   const REAL* A,
                   std::ptrdiff_t lda,
                   const REAL* x,
                   const REAL* b,
                   Rmidrad<REAL>* out);

template <class REAL>
REAL Rsum(std::ptrdiff_t n, const REAL* x, std::ptrdiff_t incx) {
    using A = Rarith<REAL>;

    if (n <= 0) {
        return A::zero();
    }
    assert(x != nullptr);
    assert(incx >= 1);

    REAL s = x[0];
    REAL c = A::zero();
    std::ptrdiff_t ix = incx;
    for (std::ptrdiff_t i = 1; i < n; ++i) {
        REAL e = A::zero();
        Rtwosum(s, x[ix], s, e);
        REAL cn = c + e;
        c = cn;
        ix += incx;
    }
    REAL result = s + c;
    return result;
}

template <class REAL>
REAL Rdot(std::ptrdiff_t n,
          const REAL* x,
          std::ptrdiff_t incx,
          const REAL* y,
          std::ptrdiff_t incy) {
    using A = Rarith<REAL>;

    if (n <= 0) {
        return A::zero();
    }
    assert(x != nullptr);
    assert(y != nullptr);
    assert(incx >= 1);
    assert(incy >= 1);

    REAL p = A::zero();
    REAL s = A::zero();
    Rtwoproduct(x[0], y[0], p, s);

    std::ptrdiff_t ix = incx;
    std::ptrdiff_t iy = incy;
    for (std::ptrdiff_t i = 1; i < n; ++i) {
        REAL h = A::zero();
        REAL r = A::zero();
        Rtwoproduct(x[ix], y[iy], h, r);
        REAL q = A::zero();
        Rtwosum(p, h, p, q);
        REAL t = q + r;
        REAL sn = s + t;
        s = sn;
        ix += incx;
        iy += incy;
    }

    REAL result = p + s;
    return result;
}

namespace detail {

template <class REAL>
Rmidrad<REAL> invalid_midrad() {
    using A = Rarith<REAL>;
    return {A::zero(), A::infinity(), Rstatus::invalid_input};
}

template <class REAL>
Rmidrad<REAL> non_finite_midrad() {
    using A = Rarith<REAL>;
    return {A::zero(), A::infinity(), Rstatus::non_finite};
}

inline Rstatus worst_status(Rstatus a, Rstatus b) {
    return (Rstatus_rank(a) >= Rstatus_rank(b)) ? a : b;
}

} // namespace detail

template <class REAL>
Rmidrad<REAL> vRsum(std::ptrdiff_t n, const REAL* x, std::ptrdiff_t incx) {
    using A = Rarith<REAL>;

    if (n < 0) {
        return detail::invalid_midrad<REAL>();
    }
    if (n == 0) {
        return {A::zero(), A::zero(), Rstatus::ok};
    }
    if (x == nullptr || incx < 1) {
        return detail::invalid_midrad<REAL>();
    }

    for (std::ptrdiff_t i = 0; i < n; ++i) {
        if (!A::is_finite(x[i * incx])) {
            return detail::non_finite_midrad<REAL>();
        }
    }

    REAL sup = A::zero();
    {
        typename A::round_up scope;
        REAL acc = A::zero();
        for (std::ptrdiff_t i = 0; i < n; ++i) {
            REAL next = acc + x[i * incx];
            acc = next;
        }
        sup = acc;
    }

    REAL inf = A::zero();
    {
        typename A::round_down scope;
        REAL acc = A::zero();
        for (std::ptrdiff_t i = 0; i < n; ++i) {
            REAL next = acc + x[i * incx];
            acc = next;
        }
        inf = acc;
    }

    REAL mid = Rmidpoint(inf, sup);
    return Rmake_midrad(inf, sup, mid);
}

template <class REAL>
Rmidrad<REAL> vRdot(std::ptrdiff_t n,
                    const REAL* x,
                    std::ptrdiff_t incx,
                    const REAL* y,
                    std::ptrdiff_t incy) {
    using A = Rarith<REAL>;

    if (n < 0) {
        return detail::invalid_midrad<REAL>();
    }
    if (n == 0) {
        return {A::zero(), A::zero(), Rstatus::ok};
    }
    if (x == nullptr || y == nullptr || incx < 1 || incy < 1) {
        return detail::invalid_midrad<REAL>();
    }

    for (std::ptrdiff_t i = 0; i < n; ++i) {
        if (!A::is_finite(x[i * incx]) || !A::is_finite(y[i * incy])) {
            return detail::non_finite_midrad<REAL>();
        }
    }

    REAL sup = A::zero();
    {
        typename A::round_up scope;
        REAL acc = A::zero();
        for (std::ptrdiff_t i = 0; i < n; ++i) {
            REAL prod = x[i * incx] * y[i * incy];
            REAL next = acc + prod;
            acc = next;
        }
        sup = acc;
    }

    REAL inf = A::zero();
    {
        typename A::round_down scope;
        REAL acc = A::zero();
        for (std::ptrdiff_t i = 0; i < n; ++i) {
            REAL prod = x[i * incx] * y[i * incy];
            REAL next = acc + prod;
            acc = next;
        }
        inf = acc;
    }

    REAL mid = Rmidpoint(inf, sup);
    return Rmake_midrad(inf, sup, mid);
}


template <class REAL>
Rstatus vRresidual(std::ptrdiff_t m,
                   std::ptrdiff_t n,
                   const REAL* A_data,
                   std::ptrdiff_t lda,
                   const REAL* x,
                   const REAL* b,
                   Rmidrad<REAL>* out) {
    using Arith = Rarith<REAL>;

    if (m < 0 || n < 0) {
        return Rstatus::invalid_input;
    }
    if (m == 0) {
        return Rstatus::ok;
    }
    if (out == nullptr || b == nullptr) {
        return Rstatus::invalid_input;
    }
    if (n == 0) {
        Rstatus worst = Rstatus::ok;
        for (std::ptrdiff_t i = 0; i < m; ++i) {
            if (!Arith::is_finite(b[i])) {
                out[i] = detail::non_finite_midrad<REAL>();
            } else {
                out[i] = {b[i], Arith::zero(), Rstatus::ok};
            }
            worst = detail::worst_status(worst, out[i].status);
        }
        return worst;
    }
    if (A_data == nullptr || x == nullptr || lda < n) {
        return Rstatus::invalid_input;
    }

    Rstatus worst = Rstatus::ok;
    for (std::ptrdiff_t row = 0; row < m; ++row) {
        bool finite = Arith::is_finite(b[row]);
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            if (!Arith::is_finite(A_data[row * lda + col]) || !Arith::is_finite(x[col])) {
                finite = false;
            }
        }
        if (!finite) {
            out[row] = detail::non_finite_midrad<REAL>();
            worst = detail::worst_status(worst, out[row].status);
            continue;
        }

        REAL sup = Arith::zero();
        {
            typename Arith::round_up scope;
            REAL acc = b[row];
            for (std::ptrdiff_t col = 0; col < n; ++col) {
                REAL na = -A_data[row * lda + col];
                REAL prod = na * x[col];
                REAL next = acc + prod;
                acc = next;
            }
            sup = acc;
        }

        REAL inf = Arith::zero();
        {
            typename Arith::round_down scope;
            REAL acc = b[row];
            for (std::ptrdiff_t col = 0; col < n; ++col) {
                REAL na = -A_data[row * lda + col];
                REAL prod = na * x[col];
                REAL next = acc + prod;
                acc = next;
            }
            inf = acc;
        }

        REAL mid = Rmidpoint(inf, sup);
        out[row] = Rmake_midrad(inf, sup, mid);
        worst = detail::worst_status(worst, out[row].status);
    }

    return worst;
}

extern template float Rsum<float>(std::ptrdiff_t, const float*, std::ptrdiff_t);
extern template double Rsum<double>(std::ptrdiff_t, const double*, std::ptrdiff_t);
extern template float Rdot<float>(std::ptrdiff_t, const float*, std::ptrdiff_t, const float*, std::ptrdiff_t);
extern template double Rdot<double>(std::ptrdiff_t, const double*, std::ptrdiff_t, const double*, std::ptrdiff_t);
extern template Rmidrad<float> vRsum<float>(std::ptrdiff_t, const float*, std::ptrdiff_t);
extern template Rmidrad<double> vRsum<double>(std::ptrdiff_t, const double*, std::ptrdiff_t);
extern template Rmidrad<float> vRdot<float>(std::ptrdiff_t, const float*, std::ptrdiff_t, const float*, std::ptrdiff_t);
extern template Rmidrad<double> vRdot<double>(std::ptrdiff_t, const double*, std::ptrdiff_t, const double*, std::ptrdiff_t);
extern template Rstatus vRresidual<float>(std::ptrdiff_t,
                                          std::ptrdiff_t,
                                          const float*,
                                          std::ptrdiff_t,
                                          const float*,
                                          const float*,
                                          Rmidrad<float>*);
extern template Rstatus vRresidual<double>(std::ptrdiff_t,
                                           std::ptrdiff_t,
                                           const double*,
                                           std::ptrdiff_t,
                                           const double*,
                                           const double*,
                                           Rmidrad<double>*);

#ifdef VMPLAPACK_ENABLE_MPFR
extern template mpfrxx::mpfr_class Rsum<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                            const mpfrxx::mpfr_class*,
                                                            std::ptrdiff_t);
extern template mpfrxx::mpfr_class Rdot<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                            const mpfrxx::mpfr_class*,
                                                            std::ptrdiff_t,
                                                            const mpfrxx::mpfr_class*,
                                                            std::ptrdiff_t);
extern template Rmidrad<mpfrxx::mpfr_class> vRsum<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                                      const mpfrxx::mpfr_class*,
                                                                      std::ptrdiff_t);
extern template Rmidrad<mpfrxx::mpfr_class> vRdot<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                                      const mpfrxx::mpfr_class*,
                                                                      std::ptrdiff_t,
                                                                      const mpfrxx::mpfr_class*,
                                                                      std::ptrdiff_t);
extern template Rstatus vRresidual<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                       std::ptrdiff_t,
                                                       const mpfrxx::mpfr_class*,
                                                       std::ptrdiff_t,
                                                       const mpfrxx::mpfr_class*,
                                                       const mpfrxx::mpfr_class*,
                                                       Rmidrad<mpfrxx::mpfr_class>*);
#endif

} // namespace vmplapack
