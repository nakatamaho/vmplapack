// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

namespace {

template <class REAL>
int output_digits() {
    return std::numeric_limits<REAL>::max_digits10;
}

template <>
int output_digits<mpfrxx::mpfr_class>() {
    return 90;
}

const char* status_name(vmplapack::Rstatus status) {
    switch (status) {
        case vmplapack::Rstatus::ok:
            return "ok";
        case vmplapack::Rstatus::unbounded:
            return "unbounded";
        case vmplapack::Rstatus::non_finite:
            return "non_finite";
        case vmplapack::Rstatus::invalid_input:
            return "invalid_input";
    }
    return "invalid_input";
}

const char* status_name(vmplapack::VerificationStatus status) {
    switch (status) {
        case vmplapack::VerificationStatus::Verified:
            return "Verified";
        case vmplapack::VerificationStatus::Unverified:
            return "Unverified";
        case vmplapack::VerificationStatus::InvalidInput:
            return "InvalidInput";
        case vmplapack::VerificationStatus::Unsupported:
            return "Unsupported";
    }
    return "Unsupported";
}

template <class REAL>
REAL power_of_two_value(int exponent) {
    if constexpr (std::is_same<REAL, float>::value || std::is_same<REAL, double>::value) {
        REAL value = static_cast<REAL>(std::ldexp(1.0, exponent));
        return value;
    } else if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        REAL value = REAL::with_precision(static_cast<mpfr_prec_t>(vmplapack::Rarith<REAL>::precision_bits()));
        mpfr_set_ui_2exp(value.mpfr_data(), 1, static_cast<mpfr_exp_t>(exponent), MPFR_RNDN);
        return value;
    } else {
        static_assert(vmplapack::Ralways_false<REAL>::value, "power_of_two_value does not support this REAL type.");
    }
}

template <class REAL>
REAL lower_display(const vmplapack::Rmidrad<REAL>& box) {
    typename vmplapack::Rarith<REAL>::round_down scope;
    REAL lo = box.mid - box.rad;
    return lo;
}

template <class REAL>
REAL upper_display(const vmplapack::Rmidrad<REAL>& box) {
    typename vmplapack::Rarith<REAL>::round_up scope;
    REAL hi = box.mid + box.rad;
    return hi;
}

template <class REAL>
void show_box(const vmplapack::Rmidrad<REAL>& box,
              const REAL& expected,
              std::ptrdiff_t row,
              std::ptrdiff_t col) {
    REAL lo = lower_display(box);
    REAL hi = upper_display(box);
    bool covered = (box.status == vmplapack::Rstatus::ok) && (lo <= expected) && (expected <= hi);

    std::cout << "  C(" << row << "," << col << ")" << '\n';
    std::cout << "    expected exact = " << expected << '\n';
    std::cout << "    mid = " << box.mid << '\n';
    std::cout << "    rad = " << box.rad << '\n';
    std::cout << "    lower display = " << lo << '\n';
    std::cout << "    upper display = " << hi << '\n';
    std::cout << "    component status = " << status_name(box.status) << '\n';
    std::cout << "    expected covered by display = " << std::boolalpha << covered << '\n';
}

template <class REAL>
void show_tier(const char* tier, int exponent) {
    using A = vmplapack::Rarith<REAL>;

    REAL one = A::one();
    REAL two = one + one;
    REAL four = two + two;
    REAL big = power_of_two_value<REAL>(exponent);

    std::vector<REAL> left;
    left.push_back(big);
    left.push_back(one);
    left.push_back(-big);
    left.push_back(one);
    left.push_back(-two);
    left.push_back(one);

    std::vector<REAL> right;
    right.push_back(one);
    right.push_back(one);
    right.push_back(one);
    right.push_back(-one);
    right.push_back(one);
    right.push_back(one);

    std::vector<REAL> expected;
    expected.push_back(one);
    expected.push_back(-one);
    expected.push_back(A::zero());
    expected.push_back(four);

    std::vector<vmplapack::Rmidrad<REAL>> result(static_cast<std::size_t>(4));
    vmplapack::VerificationStatus status =
        vmplapack::vRgemm_point<REAL>(2, 2, 3, left.data(), 3, right.data(), 2, result.data(), 2);

    std::cout << std::setprecision(output_digits<REAL>()) << std::boolalpha;
    std::cout << tier << '\n';
    std::cout << "  A = [[2^" << exponent << ", 1, -2^" << exponent << "], [1, -2, 1]]" << '\n';
    std::cout << "  B = [[1, 1], [1, -1], [1, 1]]" << '\n';
    std::cout << "  exact C = [[1, -1], [0, 4]]" << '\n';
    std::cout << "  return status = " << status_name(status) << '\n';
    std::cout << "  return status_rank = " << vmplapack::VerificationStatus_rank(status) << '\n';
    show_box(result[0], expected[0], 0, 0);
    show_box(result[1], expected[1], 0, 1);
    show_box(result[2], expected[2], 1, 0);
    show_box(result[3], expected[3], 1, 1);
}

} // namespace

int main() {
    show_tier<float>("float", 30);
    show_tier<double>("double", 60);
    mpfrxx::set_default_precision_bits(512);
    show_tier<mpfrxx::mpfr_class>("mpfr@512", 700);
    return EXIT_SUCCESS;
}
