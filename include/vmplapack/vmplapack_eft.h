// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <type_traits>

#include <vmplapack/vmplapack_arith.h>

namespace vmplapack {

template <class>
struct Ralways_false : std::false_type {};

template <class REAL>
REAL Rsplit_factor() {
    if constexpr (std::is_same<REAL, float>::value) {
        return 0x1p12f + 1.0f;
    } else if constexpr (std::is_same<REAL, double>::value) {
        return 0x1p27 + 1.0;
#ifdef VMPLAPACK_ENABLE_MPFR
    } else if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        long p = Rarith<REAL>::precision_bits();
        long half_bits = (p + 1L) / 2L;
        REAL factor = REAL::with_precision(static_cast<mpfr_prec_t>(p));
        mpfr_set_ui_2exp(factor.mpfr_data(), 1, static_cast<mpfr_exp_t>(half_bits), MPFR_RNDN);
        REAL one = Rarith<REAL>::one();
        REAL result = factor + one;
        return result;
#endif
    } else {
        static_assert(Ralways_false<REAL>::value, "Rsplit_factor is not specialized for this REAL type.");
    }
}

template <class REAL>
void Rsplit(REAL a, REAL& hi, REAL& lo) {
    REAL factor = Rsplit_factor<REAL>();
    REAL c = factor * a;
    REAL abig = c - a;
    hi = c - abig;
    lo = a - hi;
}

template <class REAL>
void Rtwosum(REAL a, REAL b, REAL& s, REAL& e) {
    s = a + b;
    REAL z = s - a;
    REAL t1 = s - z;
    REAL t2 = a - t1;
    REAL t3 = b - z;
    e = t2 + t3;
}

template <class REAL>
void Rfasttwosum(REAL a, REAL b, REAL& s, REAL& e) {
    s = a + b;
    REAL z = s - a;
    e = b - z;
}

template <class REAL>
void Rtwoproduct(REAL a, REAL b, REAL& p, REAL& e) {
    p = a * b;
#if defined(VMPLAPACK_USE_DEKKER_TWOPRODUCT)
    REAL a_hi;
    REAL a_lo;
    REAL b_hi;
    REAL b_lo;
    Rsplit(a, a_hi, a_lo);
    Rsplit(b, b_hi, b_lo);

    REAL a_hi_b_hi = a_hi * b_hi;
    REAL err1 = p - a_hi_b_hi;
    REAL a_lo_b_hi = a_lo * b_hi;
    REAL err2 = err1 - a_lo_b_hi;
    REAL a_hi_b_lo = a_hi * b_lo;
    REAL err3 = err2 - a_hi_b_lo;
    REAL a_lo_b_lo = a_lo * b_lo;
    e = a_lo_b_lo - err3;
#else
    REAL neg_p = -p;
    e = Rarith<REAL>::fma(a, b, neg_p);
#endif
}

} // namespace vmplapack
