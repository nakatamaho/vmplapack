// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#ifndef VMPLAPACK_ENABLE_MPFR
#error "Rgendot.h requires VMPLAPACK_ENABLE_MPFR because M3 generator tests use MPFR exact checks."
#endif

#include <vmplapack/vmplapack.h>

#include <cmath>
#include <type_traits>
#include <vector>

namespace vmplapack {
namespace gendot {

template <class REAL>
struct Rdot_case {
    std::vector<REAL> x;
    std::vector<REAL> y;
    REAL exact;
};

template <class REAL>
REAL power_of_two(int exponent) {
    if constexpr (std::is_same<REAL, float>::value || std::is_same<REAL, double>::value) {
        REAL one = Rarith<REAL>::one();
        REAL value = std::ldexp(one, exponent);
        return value;
    } else if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        long precision = Rarith<REAL>::precision_bits();
        REAL value = REAL::with_precision(static_cast<mpfr_prec_t>(precision));
        mpfr_set_ui_2exp(value.mpfr_data(), 1, static_cast<mpfr_exp_t>(exponent), MPFR_RNDN);
        return value;
    } else {
        static_assert(Ralways_false<REAL>::value, "gendot::power_of_two does not support this REAL type.");
    }
}

template <class REAL>
Rdot_case<REAL> family_a_alternating(int pairs, REAL delta) {
    Rdot_case<REAL> result;
    result.x.reserve(static_cast<std::size_t>(2 * pairs + 1));
    result.y.reserve(static_cast<std::size_t>(2 * pairs + 1));
    for (int i = 0; i < pairs; ++i) {
        result.x.push_back(Rarith<REAL>::one());
        result.y.push_back(Rarith<REAL>::one());
        result.x.push_back(Rarith<REAL>::one());
        result.y.push_back(-Rarith<REAL>::one());
    }
    result.x.push_back(Rarith<REAL>::one());
    result.y.push_back(delta);
    result.exact = delta;
    return result;
}

template <class REAL>
Rdot_case<REAL> family_b_two_term(REAL y2) {
    Rdot_case<REAL> result;
    result.x.push_back(Rarith<REAL>::one());
    result.x.push_back(Rarith<REAL>::one());
    result.y.push_back(Rarith<REAL>::one());
    result.y.push_back(y2);
    REAL exact = Rarith<REAL>::one() + y2;
    result.exact = exact;
    return result;
}

template <class REAL>
Rdot_case<REAL> family_c_exponent_cancellation(int exponent, REAL s) {
    Rdot_case<REAL> result;
    REAL large = power_of_two<REAL>(exponent);
    result.x.push_back(large);
    result.y.push_back(Rarith<REAL>::one());
    result.x.push_back(s);
    result.y.push_back(Rarith<REAL>::one());
    result.x.push_back(large);
    result.y.push_back(-Rarith<REAL>::one());
    result.exact = s;
    return result;
}

#ifdef VMPLAPACK_ENABLE_MPFR
inline bool has_precision(const mpfrxx::mpfr_class& value, long precision) {
    return mpfr_get_prec(value.mpfr_data()) == static_cast<mpfr_prec_t>(precision);
}

inline bool all_inputs_have_precision(const Rdot_case<mpfrxx::mpfr_class>& c, long precision) {
    if (!has_precision(c.exact, precision)) {
        return false;
    }
    for (std::size_t i = 0; i < c.x.size(); ++i) {
        if (!has_precision(c.x[i], precision) || !has_precision(c.y[i], precision)) {
            return false;
        }
    }
    return true;
}
#endif

} // namespace gendot
} // namespace vmplapack
