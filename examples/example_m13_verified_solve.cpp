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
    return 30;
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
REAL from_int(int value) {
    REAL result = static_cast<REAL>(value);
    return result;
}

template <>
mpfrxx::mpfr_class from_int<mpfrxx::mpfr_class>(int value) {
    mpfrxx::mpfr_class result =
        mpfrxx::mpfr_class::with_precision(static_cast<mpfr_prec_t>(vmplapack::Rarith<mpfrxx::mpfr_class>::precision_bits()));
    mpfr_set_si(result.mpfr_data(), static_cast<long>(value), MPFR_RNDN);
    return result;
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
void print_vector(const std::vector<REAL>& values) {
    std::cout << "[ ";
    for (std::size_t i = 0; i < values.size(); ++i) {
        std::cout << values[i];
        if (i + 1U < values.size()) {
            std::cout << ", ";
        }
    }
    std::cout << " ]";
}

template <class REAL>
void print_matrix(std::ptrdiff_t rows, std::ptrdiff_t cols, const std::vector<REAL>& values, std::ptrdiff_t ld) {
    std::cout << "[ ";
    for (std::ptrdiff_t row = 0; row < rows; ++row) {
        std::cout << "[ ";
        for (std::ptrdiff_t col = 0; col < cols; ++col) {
            std::cout << values[static_cast<std::size_t>(row * ld + col)];
            if (col + 1 < cols) {
                std::cout << ", ";
            }
        }
        if (row + 1 < rows) {
            std::cout << " ]; ";
        } else {
            std::cout << " ] ";
        }
    }
    std::cout << "]";
}

template <class REAL>
void collect_mid_radius(std::ptrdiff_t rows,
                        std::ptrdiff_t cols,
                        const std::vector<vmplapack::Rmidrad<REAL>>& boxes,
                        std::ptrdiff_t ld,
                        std::vector<REAL>& mid,
                        std::vector<REAL>& rad) {
    using A = vmplapack::Rarith<REAL>;
    mid.assign(static_cast<std::size_t>(rows * cols), A::zero());
    rad.assign(static_cast<std::size_t>(rows * cols), A::zero());
    for (std::ptrdiff_t row = 0; row < rows; ++row) {
        for (std::ptrdiff_t col = 0; col < cols; ++col) {
            const vmplapack::Rmidrad<REAL>& box = boxes[static_cast<std::size_t>(row * ld + col)];
            mid[static_cast<std::size_t>(row * cols + col)] = box.mid;
            rad[static_cast<std::size_t>(row * cols + col)] = box.rad;
        }
    }
}

template <class REAL>
REAL max_radius(std::ptrdiff_t count, const std::vector<vmplapack::Rmidrad<REAL>>& boxes, std::ptrdiff_t inc) {
    using A = vmplapack::Rarith<REAL>;
    REAL value = A::zero();
    for (std::ptrdiff_t i = 0; i < count; ++i) {
        const vmplapack::Rmidrad<REAL>& box = boxes[static_cast<std::size_t>(i * inc)];
        if (box.rad > value) {
            value = box.rad;
        }
    }
    return value;
}

template <class REAL>
REAL max_matrix_radius(std::ptrdiff_t rows,
                       std::ptrdiff_t cols,
                       const std::vector<vmplapack::Rmidrad<REAL>>& boxes,
                       std::ptrdiff_t ld) {
    using A = vmplapack::Rarith<REAL>;
    REAL value = A::zero();
    for (std::ptrdiff_t row = 0; row < rows; ++row) {
        for (std::ptrdiff_t col = 0; col < cols; ++col) {
            const vmplapack::Rmidrad<REAL>& box = boxes[static_cast<std::size_t>(row * ld + col)];
            if (box.rad > value) {
                value = box.rad;
            }
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
bool box_covers_value(const vmplapack::Rmidrad<REAL>& box, const REAL& expected) {
    if (box.status != vmplapack::Rstatus::ok) {
        return false;
    }
    REAL lo = lower_display(box);
    REAL hi = upper_display(box);
    return lo <= expected && expected <= hi;
}

template <class REAL>
void print_solution_vector(const std::vector<vmplapack::Rmidrad<REAL>>& boxes,
                           std::ptrdiff_t n,
                           std::ptrdiff_t inc,
                           const std::vector<REAL>& expected) {
    for (std::ptrdiff_t i = 0; i < n; ++i) {
        const vmplapack::Rmidrad<REAL>& box = boxes[static_cast<std::size_t>(i * inc)];
        REAL lo = lower_display(box);
        REAL hi = upper_display(box);
        std::cout << "    x(" << i << ") mid=" << box.mid
                  << " rad=" << box.rad
                  << " interval=[" << lo << ", " << hi << "]"
                  << " status=" << status_name(box.status)
                  << " expected_covered=" << std::boolalpha << box_covers_value(box, expected[static_cast<std::size_t>(i)])
                  << '\n';
    }
}

template <class REAL>
void print_matrix_coverage(std::ptrdiff_t rows,
                           std::ptrdiff_t cols,
                           const std::vector<vmplapack::Rmidrad<REAL>>& boxes,
                           std::ptrdiff_t ld,
                           const std::vector<REAL>& expected,
                           std::ptrdiff_t expected_ld) {
    for (std::ptrdiff_t row = 0; row < rows; ++row) {
        for (std::ptrdiff_t col = 0; col < cols; ++col) {
            const vmplapack::Rmidrad<REAL>& box = boxes[static_cast<std::size_t>(row * ld + col)];
            const REAL& exact = expected[static_cast<std::size_t>(row * expected_ld + col)];
            std::cout << "    X(" << row << "," << col << ") status=" << status_name(box.status)
                      << " expected_covered=" << std::boolalpha << box_covers_value(box, exact) << '\n';
        }
    }
}

template <class REAL>
void show_small_vector_solve() {
    std::vector<REAL> matrix;
    matrix.push_back(from_int<REAL>(2));
    matrix.push_back(from_int<REAL>(1));
    matrix.push_back(from_int<REAL>(-1));
    matrix.push_back(from_int<REAL>(4));

    std::vector<REAL> rhs;
    rhs.push_back(from_int<REAL>(10));
    rhs.push_back(from_int<REAL>(13));

    std::vector<REAL> expected;
    expected.push_back(from_int<REAL>(3));
    expected.push_back(from_int<REAL>(4));

    std::vector<vmplapack::Rmidrad<REAL>> result(static_cast<std::size_t>(2));
    vmplapack::VerificationStatus status = vmplapack::vRgesv<REAL>(2, matrix.data(), 2, rhs.data(), 1, result.data(), 1);
    REAL radius = max_radius(2, result, 1);

    std::cout << "  1. Small exact vector solve" << '\n';
    std::cout << "    A = ";
    print_matrix(2, 2, matrix, 2);
    std::cout << '\n';
    std::cout << "    b = ";
    print_vector(rhs);
    std::cout << '\n';
    std::cout << "    exact x = ";
    print_vector(expected);
    std::cout << '\n';
    std::cout << "    return status = " << status_name(status) << '\n';
    std::cout << "    max radius = " << radius << '\n';
    std::cout << "    max radius/u = " << radius_units(radius) << '\n';
    print_solution_vector(result, 2, 1, expected);
}

template <class REAL>
void show_matrix_rhs_solve() {
    std::vector<REAL> matrix;
    matrix.push_back(from_int<REAL>(2));
    matrix.push_back(from_int<REAL>(1));
    matrix.push_back(from_int<REAL>(-1));
    matrix.push_back(from_int<REAL>(4));

    std::vector<REAL> rhs;
    rhs.push_back(from_int<REAL>(10));
    rhs.push_back(from_int<REAL>(1));
    rhs.push_back(from_int<REAL>(13));
    rhs.push_back(from_int<REAL>(22));

    std::vector<REAL> expected;
    expected.push_back(from_int<REAL>(3));
    expected.push_back(from_int<REAL>(-2));
    expected.push_back(from_int<REAL>(4));
    expected.push_back(from_int<REAL>(5));

    std::vector<vmplapack::Rmidrad<REAL>> result(static_cast<std::size_t>(4));
    vmplapack::VerificationStatus status = vmplapack::vRgesv<REAL>(2, 2, matrix.data(), 2, rhs.data(), 2, result.data(), 2);
    std::vector<REAL> mid;
    std::vector<REAL> rad;
    collect_mid_radius(2, 2, result, 2, mid, rad);
    REAL radius = max_matrix_radius(2, 2, result, 2);

    std::cout << "  2. Matrix RHS solve" << '\n';
    std::cout << "    B = ";
    print_matrix(2, 2, rhs, 2);
    std::cout << '\n';
    std::cout << "    exact X = ";
    print_matrix(2, 2, expected, 2);
    std::cout << '\n';
    std::cout << "    return status = " << status_name(status) << '\n';
    std::cout << "    X.mid = ";
    print_matrix(2, 2, mid, 2);
    std::cout << '\n';
    std::cout << "    X.rad = ";
    print_matrix(2, 2, rad, 2);
    std::cout << '\n';
    std::cout << "    max radius/u = " << radius_units(radius) << '\n';
    print_matrix_coverage(2, 2, result, 2, expected, 2);
}

template <class REAL>
void show_near_singular_verified(int exponent) {
    using A = vmplapack::Rarith<REAL>;
    REAL one = A::one();
    REAL eps = power_of_two_value<REAL>(exponent);
    REAL one_plus_eps = one + eps;
    REAL minus_eps = -eps;

    std::vector<REAL> matrix;
    matrix.push_back(one);
    matrix.push_back(one);
    matrix.push_back(one);
    matrix.push_back(one_plus_eps);

    std::vector<REAL> rhs;
    rhs.push_back(A::zero());
    rhs.push_back(minus_eps);

    std::vector<REAL> expected;
    expected.push_back(one);
    expected.push_back(-one);

    std::vector<vmplapack::Rmidrad<REAL>> result(static_cast<std::size_t>(2));
    vmplapack::VerificationStatus status = vmplapack::vRgesv<REAL>(2, matrix.data(), 2, rhs.data(), 1, result.data(), 1);
    REAL radius = max_radius(2, result, 1);

    std::cout << "  3. Near-singular but certified" << '\n';
    std::cout << "    A = [[1, 1], [1, 1 + 2^" << exponent << "]]" << '\n';
    std::cout << "    b = [0, -2^" << exponent << "]" << '\n';
    std::cout << "    exact x = [1, -1]" << '\n';
    std::cout << "    gap = " << eps << '\n';
    std::cout << "    rough 1/gap = " << (one / eps) << '\n';
    std::cout << "    return status = " << status_name(status) << '\n';
    std::cout << "    max radius = " << radius << '\n';
    std::cout << "    max radius/u = " << radius_units(radius) << '\n';
    print_solution_vector(result, 2, 1, expected);
}

template <class REAL>
void show_residual_display() {
    std::vector<REAL> matrix;
    matrix.push_back(from_int<REAL>(2));
    matrix.push_back(from_int<REAL>(1));
    matrix.push_back(from_int<REAL>(-1));
    matrix.push_back(from_int<REAL>(4));

    std::vector<REAL> rhs;
    rhs.push_back(from_int<REAL>(10));
    rhs.push_back(from_int<REAL>(13));

    std::vector<vmplapack::Rmidrad<REAL>> solution(static_cast<std::size_t>(2));
    vmplapack::VerificationStatus solve_status = vmplapack::vRgesv<REAL>(2, matrix.data(), 2, rhs.data(), 1, solution.data(), 1);
    std::vector<REAL> midpoint;
    midpoint.push_back(solution[0].mid);
    midpoint.push_back(solution[1].mid);

    std::vector<vmplapack::Rmidrad<REAL>> residual(static_cast<std::size_t>(2));
    vmplapack::Rstatus residual_status = vmplapack::vRresidual<REAL>(2, 2, matrix.data(), 2, midpoint.data(), rhs.data(), residual.data());

    std::cout << "  4. Residual of the returned midpoint" << '\n';
    std::cout << "    solve status = " << status_name(solve_status) << '\n';
    std::cout << "    x.mid = ";
    print_vector(midpoint);
    std::cout << '\n';
    std::cout << "    residual status = " << status_name(residual_status) << '\n';
    for (std::ptrdiff_t row = 0; row < 2; ++row) {
        const vmplapack::Rmidrad<REAL>& box = residual[static_cast<std::size_t>(row)];
        REAL lo = lower_display(box);
        REAL hi = upper_display(box);
        std::cout << "    residual(" << row << ") mid=" << box.mid
                  << " rad=" << box.rad
                  << " interval=[" << lo << ", " << hi << "]"
                  << " status=" << status_name(box.status) << '\n';
    }
}

template <class REAL>
void show_certificate_failure() {
    std::vector<REAL> matrix;
    matrix.push_back(from_int<REAL>(1));
    matrix.push_back(from_int<REAL>(2));
    matrix.push_back(from_int<REAL>(2));
    matrix.push_back(from_int<REAL>(4));

    std::vector<REAL> rhs;
    rhs.push_back(from_int<REAL>(1));
    rhs.push_back(from_int<REAL>(2));

    std::vector<vmplapack::Rmidrad<REAL>> result(static_cast<std::size_t>(2));
    vmplapack::VerificationStatus status = vmplapack::vRgesv<REAL>(2, matrix.data(), 2, rhs.data(), 1, result.data(), 1);

    std::cout << "  5. Certificate failure" << '\n';
    std::cout << "    A = ";
    print_matrix(2, 2, matrix, 2);
    std::cout << '\n';
    std::cout << "    b = ";
    print_vector(rhs);
    std::cout << '\n';
    std::cout << "    return status = " << status_name(status) << '\n';
    std::cout << "    interpretation = Unverified is a failed certificate, not a singularity theorem" << '\n';
    for (std::ptrdiff_t i = 0; i < 2; ++i) {
        const vmplapack::Rmidrad<REAL>& box = result[static_cast<std::size_t>(i)];
        std::cout << "    x(" << i << ") mid=" << box.mid
                  << " rad=" << box.rad
                  << " status=" << status_name(box.status) << '\n';
    }
}

template <class REAL>
void show_tier(const char* tier, int near_singular_exponent) {
    std::cout << std::setprecision(output_digits<REAL>()) << std::boolalpha;
    std::cout << tier << '\n';
    show_small_vector_solve<REAL>();
    show_matrix_rhs_solve<REAL>();
    show_near_singular_verified<REAL>(near_singular_exponent);
    show_residual_display<REAL>();
    show_certificate_failure<REAL>();
}

} // namespace

int main() {
    show_tier<float>("float", -10);
    show_tier<double>("double", -30);
    mpfrxx::set_default_precision_bits(512);
    show_tier<mpfrxx::mpfr_class>("mpfr@512", -200);
    return 0;
}
