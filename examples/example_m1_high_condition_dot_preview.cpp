// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

namespace {

template <class REAL>
int output_digits() {
    return std::numeric_limits<REAL>::max_digits10;
}

#ifdef VMPLAPACK_ENABLE_MPFR
template <>
int output_digits<mpfrxx::mpfr_class>() {
    return 180;
}
#endif

template <class REAL>
REAL power_of_two(int exponent) {
    return std::ldexp(vmplapack::Rarith<REAL>::one(), exponent);
}

#ifdef VMPLAPACK_ENABLE_MPFR
template <>
mpfrxx::mpfr_class power_of_two<mpfrxx::mpfr_class>(int exponent) {
    using REAL = mpfrxx::mpfr_class;
    using A = vmplapack::Rarith<REAL>;

    REAL value = REAL::with_precision(static_cast<mpfr_prec_t>(A::precision_bits()));
    mpfr_set_ui_2exp(value.mpfr_data(), 1, static_cast<mpfr_exp_t>(exponent), MPFR_RNDN);
    return value;
}
#endif

template <class REAL>
REAL naive_dot(const std::vector<REAL>& x, const std::vector<REAL>& y) {
    using A = vmplapack::Rarith<REAL>;

    REAL acc = A::zero();
    for (std::size_t i = 0; i < x.size(); ++i) {
        REAL prod = x[i] * y[i];
        REAL next = acc + prod;
        acc = next;
    }
    return acc;
}

template <class REAL>
void show_high_condition_dot(const char* tier, int exponent_min, int exponent_max, int pair_count) {
    using A = vmplapack::Rarith<REAL>;

    std::uint64_t seed = static_cast<std::uint64_t>(20260222 + A::precision_bits());
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> exponent_dist(exponent_min, exponent_max);

    std::vector<REAL> x;
    std::vector<REAL> y;
    x.reserve(static_cast<std::size_t>(2 * pair_count + 1));
    y.reserve(static_cast<std::size_t>(2 * pair_count + 1));

    REAL sum_abs = A::zero();
    for (int i = 0; i < pair_count; ++i) {
        int exponent = exponent_dist(rng);
        REAL magnitude = power_of_two<REAL>(exponent);

        x.push_back(A::one());
        y.push_back(magnitude);
        x.push_back(A::one());
        y.push_back(-magnitude);

        REAL pair_abs = magnitude + magnitude;
        REAL next_sum_abs = sum_abs + pair_abs;
        sum_abs = next_sum_abs;
    }

    x.push_back(A::one());
    y.push_back(A::one());
    REAL final_sum_abs = sum_abs + A::one();
    sum_abs = final_sum_abs;

    std::vector<std::size_t> order(y.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        order[i] = i;
    }
    std::shuffle(order.begin(), order.end(), rng);

    std::vector<REAL> shuffled_x;
    std::vector<REAL> shuffled_y;
    shuffled_x.reserve(x.size());
    shuffled_y.reserve(y.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        shuffled_x.push_back(x[order[i]]);
        shuffled_y.push_back(y[order[i]]);
    }

    REAL true_dot = A::one();
    REAL naive = naive_dot(shuffled_x, shuffled_y);
    REAL error = naive - true_dot;
    REAL two = A::one() + A::one();
    REAL doubled_sum_abs = two * sum_abs;
    REAL condition_estimate = doubled_sum_abs / A::abs(true_dot);

    std::cout << std::setprecision(output_digits<REAL>());
    std::cout << tier << '\n';
    std::cout << "  seed = " << seed << '\n';
    std::cout << "  n = " << shuffled_x.size() << '\n';
    std::cout << "  known true dot = " << true_dot << '\n';
    std::cout << "  sum abs terms = " << sum_abs << '\n';
    std::cout << "  condition estimate = " << condition_estimate << '\n';
    std::cout << "  naive dot = " << naive << '\n';
    std::cout << "  naive error = " << error << '\n';
}

} // namespace

int main() {
    show_high_condition_dot<float>("float", 24, 30, 16);
    show_high_condition_dot<double>("double", 52, 60, 24);
#ifdef VMPLAPACK_ENABLE_MPFR
    mpfrxx::set_default_precision_bits(512);
    show_high_condition_dot<mpfrxx::mpfr_class>("mpfr@512", 520, 700, 16);
#endif
    return 0;
}
