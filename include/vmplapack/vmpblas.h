// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>

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
REAL Rsum_abs_dot_upper(std::ptrdiff_t n,
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
Rmidrad<REAL> vRdot_apriori(std::ptrdiff_t n,
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

template <class REAL>
REAL Rsum_abs_dot_upper(std::ptrdiff_t n,
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

    typename A::round_up scope;
    REAL acc = A::zero();
    for (std::ptrdiff_t i = 0; i < n; ++i) {
        REAL ax = A::abs(x[i * incx]);
        REAL ay = A::abs(y[i * incy]);
        REAL prod = ax * ay;
        REAL next = acc + prod;
        acc = next;
    }
    return acc;
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

inline bool apriori_length_gate(std::ptrdiff_t n, long precision_bits) {
    if (n < 0 || precision_bits <= 0) {
        return false;
    }
    long length_digits = static_cast<long>(std::numeric_limits<std::ptrdiff_t>::digits);
    if (precision_bits >= length_digits) {
        return true;
    }
    std::ptrdiff_t limit = std::ptrdiff_t{1} << static_cast<int>(precision_bits);
    return n < limit;
}

template <class REAL>
REAL length_as_real_up(std::ptrdiff_t n) {
    return static_cast<REAL>(n);
}

#ifdef VMPLAPACK_ENABLE_MPFR
template <>
inline mpfrxx::mpfr_class length_as_real_up<mpfrxx::mpfr_class>(std::ptrdiff_t n) {
    mpfrxx::mpfr_class value =
        mpfrxx::mpfr_class::with_precision(static_cast<mpfr_prec_t>(Rarith<mpfrxx::mpfr_class>::precision_bits()));
    mpfr_set_sj(value.mpfr_data(), static_cast<intmax_t>(n), MPFR_RNDU);
    return value;
}
#endif

template <class REAL>
void clear_underflow_flag() {}

template <class REAL>
bool underflow_flag_raised() {
    return false;
}

#ifdef VMPLAPACK_ENABLE_MPFR
template <>
inline void clear_underflow_flag<mpfrxx::mpfr_class>() {
    mpfr_clear_underflow();
}

template <>
inline bool underflow_flag_raised<mpfrxx::mpfr_class>() {
    return mpfr_underflow_p() != 0;
}
#endif

template <class REAL>
Rmidrad<REAL> unbounded_midrad_with_mid(REAL mid) {
    using A = Rarith<REAL>;
    REAL m = A::is_finite(mid) ? mid : A::zero();
    return {m, A::infinity(), Rstatus::unbounded};
}

template <class REAL>
REAL five() {
    using A = Rarith<REAL>;
    REAL two = A::one() + A::one();
    REAL four = two + two;
    REAL value = four + A::one();
    return value;
}

template <class REAL>
REAL Rdot_apriori_radius(std::ptrdiff_t n, REAL mid, REAL S_up) {
    using A = Rarith<REAL>;

    // ORO Prop. 5.5 Dot2 a-priori estimate, not Algorithm 5.8 Dot2Err.
    REAL u = A::unit_roundoff();
    REAL one = A::one();

    REAL n_real = A::zero();
    REAL nu = A::zero();
    {
        typename A::round_up scope;
        n_real = length_as_real_up<REAL>(n);
        REAL nu_product = n_real * u;
        nu = nu_product;
    }

    REAL gamma_den = A::zero();
    {
        typename A::round_down scope;
        REAL den = one - nu;
        gamma_den = den;
    }
    if (gamma_den <= A::zero()) {
        return A::infinity();
    }

    REAL numerator = A::zero();
    {
        typename A::round_up scope;
        REAL gamma = nu / gamma_den;
        REAL gamma_sq = gamma * gamma;
        REAL gamma_term = gamma_sq * S_up;
        REAL abs_mid = A::abs(mid);
        REAL mid_term = u * abs_mid;
        REAL underflow_n = n_real * A::eta();
        REAL underflow_term = underflow_n * five<REAL>();
        REAL partial = mid_term + gamma_term;
        numerator = partial + underflow_term;
    }

    REAL final_den = A::zero();
    {
        typename A::round_down scope;
        REAL den = one - u;
        final_den = den;
    }
    if (final_den <= A::zero()) {
        return A::infinity();
    }

    typename A::round_up scope;
    REAL rad = numerator / final_den;
    return rad;
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

    REAL mid = Rsum(n, x, incx);
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

    REAL mid = Rdot(n, x, incx, y, incy);
    return Rmake_midrad(inf, sup, mid);
}


template <class REAL>
Rmidrad<REAL> vRdot_apriori(std::ptrdiff_t n,
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

    detail::clear_underflow_flag<REAL>();
    REAL mid = Rdot(n, x, incx, y, incy);
    if (detail::underflow_flag_raised<REAL>() || !A::is_finite(mid)) {
        return detail::unbounded_midrad_with_mid(mid);
    }

    if (!detail::apriori_length_gate(n, A::precision_bits())) {
        return detail::unbounded_midrad_with_mid(mid);
    }

    detail::clear_underflow_flag<REAL>();
    REAL S_up = Rsum_abs_dot_upper(n, x, incx, y, incy);
    if (detail::underflow_flag_raised<REAL>() || !A::is_finite(S_up)) {
        return detail::unbounded_midrad_with_mid(mid);
    }

    detail::clear_underflow_flag<REAL>();
    REAL rad = detail::Rdot_apriori_radius(n, mid, S_up);
    if (detail::underflow_flag_raised<REAL>() || !A::is_finite(rad)) {
        return detail::unbounded_midrad_with_mid(mid);
    }

    return {mid, rad, Rstatus::ok};
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
extern template float Rsum_abs_dot_upper<float>(std::ptrdiff_t,
                                               const float*,
                                               std::ptrdiff_t,
                                               const float*,
                                               std::ptrdiff_t);
extern template double Rsum_abs_dot_upper<double>(std::ptrdiff_t,
                                                 const double*,
                                                 std::ptrdiff_t,
                                                 const double*,
                                                 std::ptrdiff_t);
extern template Rmidrad<float> vRsum<float>(std::ptrdiff_t, const float*, std::ptrdiff_t);
extern template Rmidrad<double> vRsum<double>(std::ptrdiff_t, const double*, std::ptrdiff_t);
extern template Rmidrad<float> vRdot<float>(std::ptrdiff_t, const float*, std::ptrdiff_t, const float*, std::ptrdiff_t);
extern template Rmidrad<double> vRdot<double>(std::ptrdiff_t, const double*, std::ptrdiff_t, const double*, std::ptrdiff_t);
extern template Rmidrad<float> vRdot_apriori<float>(std::ptrdiff_t,
                                                   const float*,
                                                   std::ptrdiff_t,
                                                   const float*,
                                                   std::ptrdiff_t);
extern template Rmidrad<double> vRdot_apriori<double>(std::ptrdiff_t,
                                                     const double*,
                                                     std::ptrdiff_t,
                                                     const double*,
                                                     std::ptrdiff_t);
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
extern template mpfrxx::mpfr_class Rsum_abs_dot_upper<mpfrxx::mpfr_class>(std::ptrdiff_t,
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
extern template Rmidrad<mpfrxx::mpfr_class> vRdot_apriori<mpfrxx::mpfr_class>(std::ptrdiff_t,
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
