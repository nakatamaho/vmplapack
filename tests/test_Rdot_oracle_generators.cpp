// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include "Rdot_oracle.h"
#include "Rgendot.h"

#include <vmplapack/vmplapack.h>

#include <cstdlib>
#include <cstddef>
#include <iostream>
#include <type_traits>

namespace {

template <class REAL>
void require(bool condition, const char* tier, const char* message) {
    if (!condition) {
        std::cerr << tier << ": " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

template <class REAL>
REAL exact_epsilon_above_one() {
    using A = vmplapack::Rarith<REAL>;

    REAL two = A::one() + A::one();
    REAL step = two * A::unit_roundoff();
    REAL next = A::one() + step;
    REAL eps = next - A::one();
    return eps;
}

template <class REAL>
void require_exact_widening(const vmplapack::gendot::Rdot_case<REAL>& c, const char* tier) {
    mpfr_prec_t precision = vmplapack::oracle::initial_precision<REAL>();

    // Catches oracle inputs that are rounded again instead of exactly widened to precision P.
    for (std::size_t i = 0; i < c.x.size(); ++i) {
        require<REAL>(vmplapack::oracle::widening_is_exact(c.x[i], precision),
                      tier,
                      "x input did not widen exactly");
        require<REAL>(vmplapack::oracle::widening_is_exact(c.y[i], precision),
                      tier,
                      "y input did not widen exactly");
    }
    require<REAL>(vmplapack::oracle::widening_is_exact(c.exact, precision),
                  tier,
                  "exact result did not widen exactly");
}

template <class REAL>
void require_oracle_stability(const vmplapack::gendot::Rdot_case<REAL>& c, const char* tier) {
    mpfr_prec_t precision = vmplapack::oracle::initial_precision<REAL>();

    // Catches an oracle precision that is too low for the cancellation in this deterministic case.
    vmplapack::oracle::Rdot_interval p_interval =
        vmplapack::oracle::dot_interval_at_precision(static_cast<std::ptrdiff_t>(c.x.size()),
                                                     c.x.data(),
                                                     1,
                                                     c.y.data(),
                                                     1,
                                                     precision);
    vmplapack::oracle::Rdot_interval two_p_interval =
        vmplapack::oracle::dot_interval_at_precision(static_cast<std::ptrdiff_t>(c.x.size()),
                                                     c.x.data(),
                                                     1,
                                                     c.y.data(),
                                                     1,
                                                     precision * 2);
    require<REAL>(vmplapack::oracle::intervals_agree_for_tier(p_interval,
                                                              two_p_interval,
                                                              vmplapack::Rarith<REAL>::precision_bits()),
                  tier,
                  "oracle interval was not stable under P to 2P");

    // Catches a directed oracle pass that fails to include the known exact result.
    require<REAL>(vmplapack::oracle::contains_exact(two_p_interval, c.exact),
                  tier,
                  "2P oracle interval does not contain the exact result");

    vmplapack::oracle::Rdot_interval self_checked =
        vmplapack::oracle::Rdot_oracle(static_cast<std::ptrdiff_t>(c.x.size()),
                                       c.x.data(),
                                       1,
                                       c.y.data(),
                                       1);
    require<REAL>(vmplapack::oracle::contains_exact(self_checked, c.exact),
                  tier,
                  "self-checked oracle interval does not contain the exact result");
}

template <class REAL>
void require_condition_exceeds_one_over_u(const vmplapack::gendot::Rdot_case<REAL>& c,
                                          const char* tier) {
    mpfr_prec_t precision = vmplapack::oracle::initial_precision<REAL>() * 2;

    // Catches Family C cases that do not actually exceed the tier-specific 1/u threshold.
    mpfrxx::mpfr_class cond = vmplapack::oracle::condition_estimate(static_cast<std::ptrdiff_t>(c.x.size()),
                                                                    c.x.data(),
                                                                    1,
                                                                    c.y.data(),
                                                                    1,
                                                                    c.exact,
                                                                    precision);
    mpfrxx::mpfr_class threshold =
        vmplapack::oracle::power_of_two(precision, vmplapack::Rarith<REAL>::precision_bits());
    require<REAL>(cond > threshold, tier, "condition estimate did not exceed 1/u");
}

template <class REAL>
void test_case(const vmplapack::gendot::Rdot_case<REAL>& c, const char* tier) {
    require_exact_widening(c, tier);
    require_oracle_stability(c, tier);
}

template <class REAL>
void test_native_tier(const char* tier, int exponent0, int exponent1) {
    using A = vmplapack::Rarith<REAL>;

    REAL eps = exact_epsilon_above_one<REAL>();
    REAL y2 = -A::one() + eps;
    test_case(vmplapack::gendot::family_a_alternating<REAL>(8, eps), tier);
    test_case(vmplapack::gendot::family_b_two_term<REAL>(y2), tier);

    vmplapack::gendot::Rdot_case<REAL> c0 =
        vmplapack::gendot::family_c_exponent_cancellation<REAL>(exponent0, A::one());
    vmplapack::gendot::Rdot_case<REAL> c1 =
        vmplapack::gendot::family_c_exponent_cancellation<REAL>(exponent1, A::one());
    test_case(c0, tier);
    test_case(c1, tier);
    require_condition_exceeds_one_over_u(c0, tier);
    require_condition_exceeds_one_over_u(c1, tier);
}

void test_mpfr_tier(long precision, bool include_w512_cases) {
    using REAL = mpfrxx::mpfr_class;
    using A = vmplapack::Rarith<REAL>;

    mpfrxx::set_default_precision_bits(static_cast<mpfr_prec_t>(precision));
    require<REAL>(A::precision_bits() == precision, "mpfr", "default precision was not set to W");

    REAL eps = exact_epsilon_above_one<REAL>();
    REAL y2 = -A::one() + eps;
    vmplapack::gendot::Rdot_case<REAL> family_a =
        vmplapack::gendot::family_a_alternating<REAL>(8, eps);
    vmplapack::gendot::Rdot_case<REAL> family_b =
        vmplapack::gendot::family_b_two_term<REAL>(y2);

    // Catches MPFR inputs accidentally built at a precision other than the fixed test W.
    require<REAL>(vmplapack::gendot::all_inputs_have_precision(family_a, precision),
                  "mpfr",
                  "Family A inputs are not at W");
    require<REAL>(vmplapack::gendot::all_inputs_have_precision(family_b, precision),
                  "mpfr",
                  "Family B inputs are not at W");
    test_case(family_a, "mpfr");
    test_case(family_b, "mpfr");

    if (include_w512_cases) {
        const int exponents[] = {520, 600, 700};
        for (int exponent : exponents) {
            vmplapack::gendot::Rdot_case<REAL> c =
                vmplapack::gendot::family_c_exponent_cancellation<REAL>(exponent, A::one());
            require<REAL>(vmplapack::gendot::all_inputs_have_precision(c, precision),
                          "mpfr",
                          "Family C inputs are not at W");
            test_case(c, "mpfr");
            require_condition_exceeds_one_over_u(c, "mpfr");
        }
    } else {
        vmplapack::gendot::Rdot_case<REAL> c =
            vmplapack::gendot::family_c_exponent_cancellation<REAL>(60, A::one());
        require<REAL>(vmplapack::gendot::all_inputs_have_precision(c, precision),
                      "mpfr",
                      "Family C inputs are not at W");
        test_case(c, "mpfr");
        require_condition_exceeds_one_over_u(c, "mpfr");
    }
}

} // namespace

int main() {
    test_native_tier<float>("float", 24, 30);
    test_native_tier<double>("double", 53, 60);
    test_mpfr_tier(53, false);
    test_mpfr_tier(512, true);
    return EXIT_SUCCESS;
}
