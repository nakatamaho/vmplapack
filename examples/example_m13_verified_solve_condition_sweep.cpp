// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

#include <cmath>
#include <cstddef>
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
    return 24;
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
bool box_covers_value(const vmplapack::Rmidrad<REAL>& box, REAL expected) {
    if (box.status != vmplapack::Rstatus::ok) {
        return false;
    }
    REAL lo = lower_display(box);
    REAL hi = upper_display(box);
    return lo <= expected && expected <= hi;
}

template <class REAL>
REAL max_radius(const std::vector<vmplapack::Rmidrad<REAL>>& boxes) {
    using A = vmplapack::Rarith<REAL>;
    REAL value = A::zero();
    for (std::size_t i = 0; i < boxes.size(); ++i) {
        if (boxes[i].rad > value) {
            value = boxes[i].rad;
        }
    }
    return value;
}

template <class REAL>
REAL radius_units(REAL radius) {
    REAL units = radius / vmplapack::Rarith<REAL>::unit_roundoff();
    return units;
}

template <class REAL>
void run_case(int exponent) {
    using A = vmplapack::Rarith<REAL>;

    REAL one = A::one();
    REAL eps = power_of_two_value<REAL>(exponent);
    REAL one_plus_eps = one + eps;
    REAL minus_eps = -eps;
    REAL minus_one = -one;

    std::vector<REAL> matrix;
    matrix.push_back(one);
    matrix.push_back(one);
    matrix.push_back(one);
    matrix.push_back(one_plus_eps);

    std::vector<REAL> rhs;
    rhs.push_back(A::zero());
    rhs.push_back(minus_eps);

    std::vector<vmplapack::Rmidrad<REAL>> result(static_cast<std::size_t>(2));
    vmplapack::VerificationStatus status = vmplapack::vRgesv<REAL>(2, matrix.data(), 2, rhs.data(), 1, result.data(), 1);
    REAL radius = max_radius(result);
    bool gap_visible = one_plus_eps > one;
    bool covered = box_covers_value(result[0], one) && box_covers_value(result[1], minus_one);

    std::cout << "    gap = 2^" << exponent
              << " visible=" << std::boolalpha << gap_visible
              << " status=" << status_name(status)
              << " max_radius=" << radius
              << " max_radius/u=" << radius_units(radius)
              << " exact_covered=" << covered
              << " component_status=[" << status_name(result[0].status) << ", " << status_name(result[1].status) << "]"
              << '\n';
}

template <class REAL>
void show_tier(const char* tier, const int* exponents, std::size_t count) {
    std::cout << std::setprecision(output_digits<REAL>()) << std::boolalpha;
    std::cout << tier << '\n';
    std::cout << "  case: A = [[1, 1], [1, 1 + gap]], b = [0, -gap], exact x = [1, -1]" << '\n';
    std::cout << "  expected: once 1 + gap rounds back to 1, the stored point matrix is singular" << '\n';
    for (std::size_t i = 0; i < count; ++i) {
        run_case<REAL>(exponents[i]);
    }
}

} // namespace

int main() {
    const int float_exponents[] = {-4, -8, -10, -16, -24, -30};
    const int double_exponents[] = {-10, -30, -40, -52, -53, -60};
    const int mpfr_exponents[] = {-50, -200, -400, -510, -512, -530};

    show_tier<float>("float", float_exponents, sizeof(float_exponents) / sizeof(float_exponents[0]));
    show_tier<double>("double", double_exponents, sizeof(double_exponents) / sizeof(double_exponents[0]));
    mpfrxx::set_default_precision_bits(512);
    show_tier<mpfrxx::mpfr_class>("mpfr@512", mpfr_exponents, sizeof(mpfr_exponents) / sizeof(mpfr_exponents[0]));
    return 0;
}
