// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#ifndef VMPLAPACK_ENABLE_MPFR
#error "Rgendot.h requires VMPLAPACK_ENABLE_MPFR because generator tests use MPFR exact checks."
#endif

#include <vmplapack/vmplapack.h>

#include "Rdot_oracle.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
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

template <class REAL>
struct Rterm {
    REAL x;
    REAL y;
};

inline mpfr_prec_t metadata_precision(long tier_precision) {
    long chosen = std::max(512L, 4L * tier_precision + 128L);
    return static_cast<mpfr_prec_t>(chosen);
}

inline mpfrxx::mpfr_class infinity_at(mpfr_prec_t precision) {
    mpfrxx::mpfr_class value = mpfrxx::mpfr_class::with_precision(precision);
    mpfr_set_inf(value.mpfr_data(), 1);
    return value;
}

template <class REAL>
const char* tier_name() {
    if constexpr (std::is_same<REAL, float>::value) {
        return "float";
    } else if constexpr (std::is_same<REAL, double>::value) {
        return "double";
    } else if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        return "mpfr";
    } else {
        static_assert(Ralways_false<REAL>::value, "tier_name does not support this REAL type.");
    }
}

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
int max_power_exponent() {
    if constexpr (std::is_same<REAL, float>::value || std::is_same<REAL, double>::value) {
        return std::numeric_limits<REAL>::max_exponent - 2;
    } else if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        mpfr_exp_t emax = mpfr_get_emax();
        if (emax > static_cast<mpfr_exp_t>(std::numeric_limits<int>::max())) {
            return std::numeric_limits<int>::max() - 2;
        }
        return static_cast<int>(emax - 2);
    } else {
        static_assert(Ralways_false<REAL>::value, "max_power_exponent does not support this REAL type.");
    }
}

template <class REAL>
int min_power_exponent() {
    if constexpr (std::is_same<REAL, float>::value || std::is_same<REAL, double>::value) {
        return std::numeric_limits<REAL>::min_exponent + 2;
    } else if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        mpfr_exp_t emin = mpfr_get_emin();
        if (emin < static_cast<mpfr_exp_t>(std::numeric_limits<int>::min())) {
            return std::numeric_limits<int>::min() + 2;
        }
        return static_cast<int>(emin + 2);
    } else {
        static_assert(Ralways_false<REAL>::value, "min_power_exponent does not support this REAL type.");
    }
}

inline int ceil_log2_size(std::size_t n) {
    if (n <= 1U) {
        return 0;
    }
    std::size_t value = n - 1U;
    int result = 0;
    while (value != 0U) {
        value >>= 1U;
        ++result;
    }
    return result;
}

template <class REAL>
Rterm<REAL> signed_power_term(int exponent, bool negative) {
    Rterm<REAL> term;
    term.x = power_of_two<REAL>(exponent);
    term.y = negative ? -Rarith<REAL>::one() : Rarith<REAL>::one();
    return term;
}

template <class REAL>
void terms_to_case(const std::vector<Rterm<REAL>>& terms, REAL exact, Rdot_case<REAL>& c) {
    c.x.clear();
    c.y.clear();
    c.x.reserve(terms.size());
    c.y.reserve(terms.size());
    for (std::size_t i = 0; i < terms.size(); ++i) {
        c.x.push_back(terms[i].x);
        c.y.push_back(terms[i].y);
    }
    c.exact = exact;
}

template <class REAL>
bool term_abs_greater(const Rterm<REAL>& a, const Rterm<REAL>& b) {
    using A = Rarith<REAL>;

    REAL prod_a = a.x * a.y;
    REAL prod_b = b.x * b.y;
    REAL abs_a = A::abs(prod_a);
    REAL abs_b = A::abs(prod_b);
    if (abs_a > abs_b) {
        return true;
    }
    if (abs_b > abs_a) {
        return false;
    }
    return a.y > b.y;
}

template <class REAL>
void apply_permutation(std::vector<Rterm<REAL>>& terms, Rpermutation permutation, std::uint64_t seed) {
    if (permutation == Rpermutation::sorted) {
        std::sort(terms.begin(), terms.end(), term_abs_greater<REAL>);
    } else if (permutation == Rpermutation::reversed) {
        std::reverse(terms.begin(), terms.end());
    } else if (permutation == Rpermutation::shuffled) {
        std::mt19937_64 engine(seed);
        std::shuffle(terms.begin(), terms.end(), engine);
    }
}

template <class REAL>
void measure_conditions(Rgenerated_dot_case<REAL>& result) {
    long tier_precision = Rarith<REAL>::precision_bits();
    vmplapack::oracle::Rdot_interval interval =
        vmplapack::oracle::Rdot_oracle(static_cast<std::ptrdiff_t>(result.data.x.size()),
                                       result.data.x.data(),
                                       1,
                                       result.data.y.data(),
                                       1);
    mpfr_prec_t precision = interval.precision + static_cast<mpfr_prec_t>(tier_precision + 128L);
    mpfrxx::mpfr_class lo = vmplapack::oracle::widen_mpfr(interval.lo, precision);
    mpfrxx::mpfr_class hi = vmplapack::oracle::widen_mpfr(interval.hi, precision);
    mpfrxx::mpfr_class sum = mpfrxx::mpfr_class::with_precision(precision);
    sum = lo + hi;
    mpfrxx::mpfr_class half = mpfrxx::mpfr_class::with_precision(precision, 0.5);
    mpfrxx::mpfr_class dot_mid = mpfrxx::mpfr_class::with_precision(precision);
    dot_mid = sum * half;
    mpfrxx::mpfr_class abs_dot = vmplapack::oracle::abs_at(dot_mid, precision);
    mpfrxx::mpfr_class s_up = vmplapack::oracle::upward_abs_term_sum(static_cast<std::ptrdiff_t>(result.data.x.size()),
                                                                     result.data.x.data(),
                                                                     1,
                                                                     result.data.y.data(),
                                                                     1,
                                                                     precision);

    if (abs_dot == vmplapack::oracle::zero_at(precision)) {
        if (s_up == vmplapack::oracle::zero_at(precision)) {
            result.measured_cond_sum = vmplapack::oracle::zero_at(precision);
            result.measured_cond_oro = vmplapack::oracle::zero_at(precision);
        } else {
            result.measured_cond_sum = infinity_at(precision);
            result.measured_cond_oro = infinity_at(precision);
        }
        return;
    }

    {
        mpfrxx::rounding_mode_scope scope(MPFR_RNDU);
        result.measured_cond_sum = s_up / abs_dot;
        mpfrxx::mpfr_class two = mpfrxx::mpfr_class::with_precision(precision, 2.0);
        result.measured_cond_oro = two * result.measured_cond_sum;
    }
}

template <class REAL>
Rgenerated_dot_case<REAL> make_generated_result(Rdot_case<REAL> c,
                                                mpfrxx::mpfr_class target_cond_oro,
                                                Rgenerator_status status,
                                                std::uint64_t seed,
                                                int scale,
                                                Rpermutation permutation) {
    Rgenerated_dot_case<REAL> result;
    result.data = c;
    result.status = status;
    result.target_cond_oro = target_cond_oro;
    mpfr_prec_t precision = metadata_precision(Rarith<REAL>::precision_bits());
    result.measured_cond_oro = vmplapack::oracle::zero_at(precision);
    result.measured_cond_sum = vmplapack::oracle::zero_at(precision);
    result.seed = seed;
    result.scale = scale;
    result.permutation = permutation;
    result.tier = tier_name<REAL>();
    if (status == Rgenerator_status::ok) {
        measure_conditions(result);
    }
    return result;
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

template <class REAL>
Rgenerated_dot_case<REAL> randomized_high_condition_power2(int target_cond_power2,
                                                           std::uint64_t seed,
                                                           Rpermutation permutation,
                                                           std::size_t pairs = 4U) {
    mpfr_prec_t precision = metadata_precision(Rarith<REAL>::precision_bits());
    mpfrxx::mpfr_class target_cond = vmplapack::oracle::power_of_two(precision,
                                                                     static_cast<mpfr_exp_t>(target_cond_power2));
    Rdot_case<REAL> empty;
    empty.exact = Rarith<REAL>::zero();

    if (target_cond_power2 < 1 || pairs == 0U) {
        return make_generated_result(empty, target_cond, Rgenerator_status::invalid_target, seed, 0, permutation);
    }

    int scale = target_cond_power2 - ceil_log2_size(4U * pairs);
    if (scale < min_power_exponent<REAL>() || scale > max_power_exponent<REAL>()) {
        return make_generated_result(empty, target_cond, Rgenerator_status::unachievable, seed, scale, permutation);
    }

    std::vector<Rterm<REAL>> terms;
    terms.reserve(2U * pairs + 1U);
    std::mt19937_64 engine(seed);
    for (std::size_t i = 0; i < pairs; ++i) {
        bool negative_first = (engine() & 1ULL) != 0ULL;
        terms.push_back(signed_power_term<REAL>(scale, negative_first));
        terms.push_back(signed_power_term<REAL>(scale, !negative_first));
    }
    terms.push_back(signed_power_term<REAL>(0, false));
    apply_permutation(terms, permutation, seed);

    Rdot_case<REAL> generated;
    terms_to_case(terms, Rarith<REAL>::one(), generated);
    return make_generated_result(generated, target_cond, Rgenerator_status::ok, seed, scale, permutation);
}

template <class REAL>
Rgenerated_dot_case<REAL> adversarial_exact_cancellation(int scale,
                                                         std::uint64_t seed,
                                                         Rpermutation permutation) {
    std::vector<Rterm<REAL>> terms;
    terms.push_back(signed_power_term<REAL>(scale, false));
    terms.push_back(signed_power_term<REAL>(scale, true));
    apply_permutation(terms, permutation, seed);

    Rdot_case<REAL> c;
    terms_to_case(terms, Rarith<REAL>::zero(), c);
    mpfr_prec_t precision = metadata_precision(Rarith<REAL>::precision_bits());
    mpfrxx::mpfr_class target = infinity_at(precision);
    return make_generated_result(c, target, Rgenerator_status::ok, seed, scale, permutation);
}

template <class REAL>
Rgenerated_dot_case<REAL> adversarial_heavy_cancellation(int scale,
                                                         std::uint64_t seed,
                                                         Rpermutation permutation) {
    std::vector<Rterm<REAL>> terms;
    terms.push_back(signed_power_term<REAL>(scale, false));
    terms.push_back(signed_power_term<REAL>(0, false));
    terms.push_back(signed_power_term<REAL>(scale, true));
    apply_permutation(terms, permutation, seed);

    Rdot_case<REAL> c;
    terms_to_case(terms, Rarith<REAL>::one(), c);
    mpfr_prec_t precision = metadata_precision(Rarith<REAL>::precision_bits());
    mpfrxx::mpfr_class target = vmplapack::oracle::zero_at(precision);
    return make_generated_result(c, target, Rgenerator_status::ok, seed, scale, permutation);
}

template <class REAL>
Rgenerated_dot_case<REAL> adversarial_alternating_signs(int pairs,
                                                        int residual_exponent,
                                                        std::uint64_t seed,
                                                        Rpermutation permutation) {
    std::vector<Rterm<REAL>> terms;
    terms.reserve(static_cast<std::size_t>(2 * pairs + 1));
    for (int i = 0; i < pairs; ++i) {
        terms.push_back(signed_power_term<REAL>(0, false));
        terms.push_back(signed_power_term<REAL>(0, true));
    }
    terms.push_back(signed_power_term<REAL>(residual_exponent, false));
    apply_permutation(terms, permutation, seed);

    Rdot_case<REAL> c;
    terms_to_case(terms, power_of_two<REAL>(residual_exponent), c);
    mpfr_prec_t precision = metadata_precision(Rarith<REAL>::precision_bits());
    mpfrxx::mpfr_class target = vmplapack::oracle::zero_at(precision);
    return make_generated_result(c, target, Rgenerator_status::ok, seed, residual_exponent, permutation);
}

template <class REAL>
Rgenerated_dot_case<REAL> adversarial_huge_small_mixing(int large_exponent,
                                                        int small_exponent,
                                                        std::uint64_t seed,
                                                        Rpermutation permutation) {
    std::vector<Rterm<REAL>> terms;
    terms.push_back(signed_power_term<REAL>(large_exponent, false));
    terms.push_back(signed_power_term<REAL>(small_exponent, false));
    terms.push_back(signed_power_term<REAL>(large_exponent, true));
    terms.push_back(signed_power_term<REAL>(small_exponent, true));
    terms.push_back(signed_power_term<REAL>(0, false));
    apply_permutation(terms, permutation, seed);

    Rdot_case<REAL> c;
    terms_to_case(terms, Rarith<REAL>::one(), c);
    mpfr_prec_t precision = metadata_precision(Rarith<REAL>::precision_bits());
    mpfrxx::mpfr_class target = vmplapack::oracle::zero_at(precision);
    return make_generated_result(c, target, Rgenerator_status::ok, seed, large_exponent, permutation);
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
