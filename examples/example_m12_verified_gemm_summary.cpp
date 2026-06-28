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
bool covers_expected(const vmplapack::Rmidrad<REAL>& box, const REAL& expected) {
    REAL lo = lower_display(box);
    REAL hi = upper_display(box);
    return box.status == vmplapack::Rstatus::ok && lo <= expected && expected <= hi;
}

template <class REAL>
void append_problem(std::vector<REAL>& left, std::vector<REAL>& right, std::vector<REAL>& expected, int exponent) {
    using A = vmplapack::Rarith<REAL>;

    REAL big = power_of_two_value<REAL>(exponent);
    REAL one = A::one();
    REAL two = one + one;
    REAL three = two + one;
    REAL four = two + two;
    REAL five = four + one;
    REAL six = three + three;
    REAL seven = six + one;
    REAL eight = four + four;
    REAL nine = eight + one;

    left.push_back(big);
    left.push_back(one);
    left.push_back(-big);
    left.push_back(two);
    left.push_back(-three);
    left.push_back(five);

    left.push_back(big);
    left.push_back(-two);
    left.push_back(-big);
    left.push_back(-one);
    left.push_back(four);
    left.push_back(-six);

    left.push_back(one);
    left.push_back(-one);
    left.push_back(one);
    left.push_back(-one);
    left.push_back(one);
    left.push_back(-one);

    left.push_back(big);
    left.push_back(A::zero());
    left.push_back(-big);
    left.push_back(seven);
    left.push_back(-eight);
    left.push_back(nine);

    right.push_back(one);
    right.push_back(one);
    right.push_back(one);
    right.push_back(one);

    right.push_back(one);
    right.push_back(-one);
    right.push_back(A::zero());
    right.push_back(two);

    right.push_back(one);
    right.push_back(one);
    right.push_back(one);
    right.push_back(one);

    right.push_back(one);
    right.push_back(-one);
    right.push_back(A::zero());
    right.push_back(-two);

    right.push_back(one);
    right.push_back(one);
    right.push_back(one);
    right.push_back(one);

    right.push_back(one);
    right.push_back(-one);
    right.push_back(A::zero());
    right.push_back(two);

    const int values[16] = {
        5,  -11, -3,  5,
        -5, 13,  4,  -10,
        0,  6,   3,  1,
        8,  -24, -8, -4
    };
    for (std::size_t i = 0; i < 16U; ++i) {
        expected.push_back(from_int<REAL>(values[i]));
    }
}

template <class REAL>
void print_octave_matrix(std::ptrdiff_t rows, std::ptrdiff_t cols, const REAL* data, std::ptrdiff_t lda) {
    std::cout << "[ ";
    for (std::ptrdiff_t row = 0; row < rows; ++row) {
        std::cout << "[ ";
        for (std::ptrdiff_t col = 0; col < cols; ++col) {
            const REAL& value = data[static_cast<std::size_t>(row * lda + col)];
            std::cout << value;
            if (col < cols - 1) {
                std::cout << ", ";
            }
        }
        if (row < rows - 1) {
            std::cout << "]; ";
        } else {
            std::cout << "] ";
        }
    }
    std::cout << "]";
}

template <class REAL>
void print_mid_matrix(std::ptrdiff_t rows,
                      std::ptrdiff_t cols,
                      const std::vector<vmplapack::Rmidrad<REAL>>& data,
                      std::ptrdiff_t ld) {
    std::cout << "[ ";
    for (std::ptrdiff_t row = 0; row < rows; ++row) {
        std::cout << "[ ";
        for (std::ptrdiff_t col = 0; col < cols; ++col) {
            const vmplapack::Rmidrad<REAL>& value = data[static_cast<std::size_t>(row * ld + col)];
            std::cout << value.mid;
            if (col < cols - 1) {
                std::cout << ", ";
            }
        }
        if (row < rows - 1) {
            std::cout << "]; ";
        } else {
            std::cout << "] ";
        }
    }
    std::cout << "]";
}

template <class REAL>
void print_rad_matrix(std::ptrdiff_t rows,
                      std::ptrdiff_t cols,
                      const std::vector<vmplapack::Rmidrad<REAL>>& data,
                      std::ptrdiff_t ld) {
    std::cout << "[ ";
    for (std::ptrdiff_t row = 0; row < rows; ++row) {
        std::cout << "[ ";
        for (std::ptrdiff_t col = 0; col < cols; ++col) {
            const vmplapack::Rmidrad<REAL>& value = data[static_cast<std::size_t>(row * ld + col)];
            std::cout << value.rad;
            if (col < cols - 1) {
                std::cout << ", ";
            }
        }
        if (row < rows - 1) {
            std::cout << "]; ";
        } else {
            std::cout << "] ";
        }
    }
    std::cout << "]";
}

template <class REAL>
void show_tier(const char* tier, int exponent) {
    using A = vmplapack::Rarith<REAL>;

    std::vector<REAL> left;
    std::vector<REAL> right;
    std::vector<REAL> expected;
    append_problem(left, right, expected, exponent);

    std::vector<vmplapack::Rmidrad<REAL>> result(static_cast<std::size_t>(16));
    vmplapack::VerificationStatus status =
        vmplapack::vRgemm_point<REAL>(4, 4, 6, left.data(), 6, right.data(), 4, result.data(), 4);

    bool all_ok = true;
    bool all_covered = true;
    REAL max_radius = A::zero();
    std::ptrdiff_t worst_index = 0;
    for (std::ptrdiff_t i = 0; i < 16; ++i) {
        const vmplapack::Rmidrad<REAL>& box = result[static_cast<std::size_t>(i)];
        if (box.status != vmplapack::Rstatus::ok) {
            all_ok = false;
        }
        if (!covers_expected(box, expected[static_cast<std::size_t>(i)])) {
            all_covered = false;
        }
        if (box.rad > max_radius) {
            max_radius = box.rad;
            worst_index = i;
        }
    }

    std::ptrdiff_t worst_row = worst_index / 4;
    std::ptrdiff_t worst_col = worst_index % 4;
    const vmplapack::Rmidrad<REAL>& worst = result[static_cast<std::size_t>(worst_index)];
    REAL worst_lo = lower_display(worst);
    REAL worst_hi = upper_display(worst);

    std::cout << std::setprecision(output_digits<REAL>()) << std::boolalpha;
    std::cout << tier << '\n';
    std::cout << "  problem = A(4x6) * B(6x4), row-major" << '\n';
    std::cout << "  cancellation scale = 2^" << exponent << '\n';
    std::cout << "  A = ";
    print_octave_matrix(4, 6, left.data(), 6);
    std::cout << '\n';
    std::cout << "  B = ";
    print_octave_matrix(6, 4, right.data(), 4);
    std::cout << '\n';
    std::cout << "  exact C = ";
    print_octave_matrix(4, 4, expected.data(), 4);
    std::cout << '\n';
    std::cout << "  result C.mid = ";
    print_mid_matrix(4, 4, result, 4);
    std::cout << '\n';
    std::cout << "  result C.rad = ";
    print_rad_matrix(4, 4, result, 4);
    std::cout << '\n';
    std::cout << "  return status = " << status_name(status) << '\n';
    std::cout << "  all components ok = " << all_ok << '\n';
    std::cout << "  all expected values covered = " << all_covered << '\n';
    std::cout << "  max radius = " << max_radius << '\n';
    std::cout << "  worst radius component = C(" << worst_row << "," << worst_col << ")" << '\n';
    std::cout << "    expected = " << expected[static_cast<std::size_t>(worst_index)] << '\n';
    std::cout << "    mid = " << worst.mid << '\n';
    std::cout << "    interval = [" << worst_lo << ", " << worst_hi << "]" << '\n';
    std::cout << "    status = " << status_name(worst.status) << '\n';
}

} // namespace

int main() {
    show_tier<float>("float", 30);
    show_tier<double>("double", 60);
    mpfrxx::set_default_precision_bits(512);
    show_tier<mpfrxx::mpfr_class>("mpfr@512", 700);
    return EXIT_SUCCESS;
}
