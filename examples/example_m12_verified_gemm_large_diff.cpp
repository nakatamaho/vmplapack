// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include "../tests/Rdot_oracle.h"

#include <vmplapack/vmplapack.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

namespace {

constexpr std::ptrdiff_t rows = 2;
constexpr std::ptrdiff_t inner = 4;
constexpr std::ptrdiff_t cols = 2;

template <class REAL>
struct large_diff_exponents;

template <>
struct large_diff_exponents<float> {
    static constexpr int p = 100;
    static constexpr int h1 = 70;
    static constexpr int h2 = 40;
};

template <>
struct large_diff_exponents<double> {
    static constexpr int p = 900;
    static constexpr int h1 = 800;
    static constexpr int h2 = 700;
};

template <>
struct large_diff_exponents<mpfrxx::mpfr_class> {
    static constexpr int p = 1800;
    static constexpr int h1 = 1200;
    static constexpr int h2 = 600;
};

template <class REAL>
int output_digits() {
    return std::numeric_limits<REAL>::max_digits10;
}

template <>
int output_digits<mpfrxx::mpfr_class>() {
    return 20;
}

int oracle_output_digits() {
    return 20;
}

int diff_output_digits() {
    return 3;
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

mpfrxx::mpfr_class interval_midpoint(const vmplapack::oracle::Rdot_interval& interval,
                                     mpfr_prec_t precision) {
    mpfrxx::mpfr_class lo = vmplapack::oracle::widen_mpfr(interval.lo, precision);
    mpfrxx::mpfr_class hi = vmplapack::oracle::widen_mpfr(interval.hi, precision);
    mpfrxx::mpfr_class sum = mpfrxx::mpfr_class::with_precision(precision);
    sum = lo + hi;
    mpfrxx::mpfr_class half = mpfrxx::mpfr_class::with_precision(precision, 0.5);
    mpfrxx::mpfr_class mid = mpfrxx::mpfr_class::with_precision(precision);
    mid = sum * half;
    return mid;
}

template <class REAL>
bool midrad_covers_oracle(const vmplapack::Rmidrad<REAL>& box,
                          const vmplapack::oracle::Rdot_interval& ref) {
    using A = vmplapack::Rarith<REAL>;

    if (box.status != vmplapack::Rstatus::ok || box.rad < A::zero()) {
        return false;
    }

    mpfr_prec_t precision = ref.precision + static_cast<mpfr_prec_t>(A::precision_bits() + 128L);
    mpfrxx::mpfr_class mid = vmplapack::oracle::widen_value(box.mid, precision);
    mpfrxx::mpfr_class rad = vmplapack::oracle::widen_value(box.rad, precision);
    mpfrxx::mpfr_class lo = mpfrxx::mpfr_class::with_precision(precision);
    mpfrxx::mpfr_class hi = mpfrxx::mpfr_class::with_precision(precision);
    {
        mpfrxx::rounding_mode_scope scope(MPFR_RNDD);
        lo = mid - rad;
    }
    {
        mpfrxx::rounding_mode_scope scope(MPFR_RNDU);
        hi = mid + rad;
    }
    mpfrxx::mpfr_class ref_lo = vmplapack::oracle::widen_mpfr(ref.lo, precision);
    mpfrxx::mpfr_class ref_hi = vmplapack::oracle::widen_mpfr(ref.hi, precision);
    return lo <= ref_lo && ref_hi <= hi;
}

template <class REAL>
mpfrxx::mpfr_class midpoint_difference(const vmplapack::Rmidrad<REAL>& box,
                                       const mpfrxx::mpfr_class& oracle_mid,
                                       mpfr_prec_t precision) {
    mpfrxx::mpfr_class mid = vmplapack::oracle::widen_value(box.mid, precision);
    mpfrxx::mpfr_class ref = vmplapack::oracle::widen_mpfr(oracle_mid, precision);
    mpfrxx::mpfr_class diff = mpfrxx::mpfr_class::with_precision(precision);
    diff = mid - ref;
    return diff;
}

template <class REAL>
mpfrxx::mpfr_class abs_midpoint_difference(const vmplapack::Rmidrad<REAL>& box,
                                           const mpfrxx::mpfr_class& oracle_mid,
                                           mpfr_prec_t precision) {
    mpfrxx::mpfr_class diff = midpoint_difference(box, oracle_mid, precision);
    mpfrxx::mpfr_class abs_diff = vmplapack::oracle::abs_at(diff, precision);
    return abs_diff;
}

template <class REAL>
void append_problem(std::vector<REAL>& left, std::vector<REAL>& right) {
    using A = vmplapack::Rarith<REAL>;

    REAL p = power_of_two_value<REAL>(large_diff_exponents<REAL>::p);
    REAL h1 = power_of_two_value<REAL>(large_diff_exponents<REAL>::h1);
    REAL h2 = power_of_two_value<REAL>(large_diff_exponents<REAL>::h2);
    REAL one = A::one();
    REAL minus_one = -one;

    left.clear();
    right.clear();
    left.reserve(static_cast<std::size_t>(rows * inner));
    right.reserve(static_cast<std::size_t>(inner * cols));

    left.push_back(p);
    left.push_back(h1);
    left.push_back(h2);
    left.push_back(-p);

    left.push_back(p);
    left.push_back(-h1);
    left.push_back(h2);
    left.push_back(-p);

    right.push_back(one);
    right.push_back(one);

    right.push_back(one);
    right.push_back(minus_one);

    right.push_back(one);
    right.push_back(minus_one);

    right.push_back(one);
    right.push_back(one);
}

template <class REAL>
void print_octave_matrix(std::ptrdiff_t matrix_rows,
                         std::ptrdiff_t matrix_cols,
                         const REAL* data,
                         std::ptrdiff_t lda) {
    std::cout << "[ ";
    for (std::ptrdiff_t row = 0; row < matrix_rows; ++row) {
        std::cout << "[ ";
        for (std::ptrdiff_t col = 0; col < matrix_cols; ++col) {
            const REAL& value = data[static_cast<std::size_t>(row * lda + col)];
            std::cout << value;
            if (col < matrix_cols - 1) {
                std::cout << ", ";
            }
        }
        if (row < matrix_rows - 1) {
            std::cout << " ]; ";
        } else {
            std::cout << " ] ";
        }
    }
    std::cout << "]";
}

void print_mpfr_matrix(std::ptrdiff_t matrix_rows,
                       std::ptrdiff_t matrix_cols,
                       const std::vector<mpfrxx::mpfr_class>& data,
                       std::ptrdiff_t ld) {
    std::cout << "[ ";
    for (std::ptrdiff_t row = 0; row < matrix_rows; ++row) {
        std::cout << "[ ";
        for (std::ptrdiff_t col = 0; col < matrix_cols; ++col) {
            const mpfrxx::mpfr_class& value = data[static_cast<std::size_t>(row * ld + col)];
            std::cout << value;
            if (col < matrix_cols - 1) {
                std::cout << ", ";
            }
        }
        if (row < matrix_rows - 1) {
            std::cout << " ]; ";
        } else {
            std::cout << " ] ";
        }
    }
    std::cout << "]";
}

template <class REAL>
void print_mid_matrix(std::ptrdiff_t matrix_rows,
                      std::ptrdiff_t matrix_cols,
                      const std::vector<vmplapack::Rmidrad<REAL>>& data,
                      std::ptrdiff_t ld) {
    std::cout << "[ ";
    for (std::ptrdiff_t row = 0; row < matrix_rows; ++row) {
        std::cout << "[ ";
        for (std::ptrdiff_t col = 0; col < matrix_cols; ++col) {
            const vmplapack::Rmidrad<REAL>& value = data[static_cast<std::size_t>(row * ld + col)];
            std::cout << value.mid;
            if (col < matrix_cols - 1) {
                std::cout << ", ";
            }
        }
        if (row < matrix_rows - 1) {
            std::cout << " ]; ";
        } else {
            std::cout << " ] ";
        }
    }
    std::cout << "]";
}

template <class REAL>
void print_rad_matrix(std::ptrdiff_t matrix_rows,
                      std::ptrdiff_t matrix_cols,
                      const std::vector<vmplapack::Rmidrad<REAL>>& data,
                      std::ptrdiff_t ld) {
    std::cout << "[ ";
    for (std::ptrdiff_t row = 0; row < matrix_rows; ++row) {
        std::cout << "[ ";
        for (std::ptrdiff_t col = 0; col < matrix_cols; ++col) {
            const vmplapack::Rmidrad<REAL>& value = data[static_cast<std::size_t>(row * ld + col)];
            std::cout << value.rad;
            if (col < matrix_cols - 1) {
                std::cout << ", ";
            }
        }
        if (row < matrix_rows - 1) {
            std::cout << " ]; ";
        } else {
            std::cout << " ] ";
        }
    }
    std::cout << "]";
}

template <class REAL>
void show_tier(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    std::vector<REAL> left;
    std::vector<REAL> right;
    append_problem(left, right);

    std::vector<vmplapack::Rmidrad<REAL>> result(static_cast<std::size_t>(rows * cols));
    vmplapack::VerificationStatus status =
        vmplapack::vRgemm_point<REAL>(rows, cols, inner, left.data(), inner, right.data(), cols, result.data(), cols);

    bool all_ok = true;
    bool all_covered = true;
    REAL max_radius = A::zero();
    mpfr_prec_t metric_precision = static_cast<mpfr_prec_t>(std::max(1024L, 4L * A::precision_bits() + 256L));
    mpfrxx::mpfr_class max_midpoint_error = vmplapack::oracle::zero_at(metric_precision);
    std::vector<mpfrxx::mpfr_class> oracle_midpoints;
    std::vector<mpfrxx::mpfr_class> midpoint_diffs;
    oracle_midpoints.reserve(static_cast<std::size_t>(rows * cols));
    midpoint_diffs.reserve(static_cast<std::size_t>(rows * cols));

    for (std::ptrdiff_t row = 0; row < rows; ++row) {
        for (std::ptrdiff_t col = 0; col < cols; ++col) {
            std::ptrdiff_t index = row * cols + col;
            const vmplapack::Rmidrad<REAL>& box = result[static_cast<std::size_t>(index)];
            vmplapack::oracle::Rdot_interval ref =
                vmplapack::oracle::Rdot_oracle(inner, left.data() + row * inner, 1, right.data() + col, cols);
            mpfr_prec_t local_precision = ref.precision + static_cast<mpfr_prec_t>(A::precision_bits() + 128L);
            mpfrxx::mpfr_class oracle_mid = interval_midpoint(ref, local_precision);
            oracle_midpoints.push_back(oracle_mid);
            midpoint_diffs.push_back(midpoint_difference(box, oracle_mid, metric_precision));
            if (box.status != vmplapack::Rstatus::ok) {
                all_ok = false;
            }
            if (!midrad_covers_oracle(box, ref)) {
                all_covered = false;
            }
            mpfrxx::mpfr_class abs_diff = abs_midpoint_difference(box, oracle_mid, metric_precision);
            if (abs_diff > max_midpoint_error) {
                max_midpoint_error = abs_diff;
            }
            if (box.rad > max_radius) {
                max_radius = box.rad;
            }
        }
    }

    std::cout << std::setprecision(output_digits<REAL>()) << std::boolalpha;
    std::cout << tier << '\n';
    std::cout << "  construction = [2^p, 2^h1, 2^h2, -2^p] dot signed columns" << '\n';
    std::cout << "  exponents p/h1/h2 = " << large_diff_exponents<REAL>::p << "/"
              << large_diff_exponents<REAL>::h1 << "/" << large_diff_exponents<REAL>::h2 << '\n';
    std::cout << "  A = ";
    print_octave_matrix(rows, inner, left.data(), inner);
    std::cout << '\n';
    std::cout << "  B = ";
    print_octave_matrix(inner, cols, right.data(), cols);
    std::cout << '\n';
    std::cout << "  oracle C midpoint = ";
    std::cout << std::setprecision(oracle_output_digits());
    print_mpfr_matrix(rows, cols, oracle_midpoints, cols);
    std::cout << '\n';
    std::cout << std::setprecision(output_digits<REAL>());
    std::cout << "  result C.mid = ";
    print_mid_matrix(rows, cols, result, cols);
    std::cout << '\n';
    std::cout << "  result C.rad = ";
    print_rad_matrix(rows, cols, result, cols);
    std::cout << '\n';
    std::cout << "  result C.diff = ";
    std::cout << std::scientific << std::setprecision(diff_output_digits());
    print_mpfr_matrix(rows, cols, midpoint_diffs, cols);
    std::cout << '\n';
    std::cout << std::defaultfloat << std::setprecision(output_digits<REAL>());
    std::cout << "  return status = " << status_name(status) << '\n';
    std::cout << "  all components ok = " << all_ok << '\n';
    std::cout << "  all oracle intervals covered = " << all_covered << '\n';
    std::cout << std::setprecision(oracle_output_digits());
    std::cout << "  max |C.mid - oracle midpoint| = " << max_midpoint_error << '\n';
    std::cout << std::setprecision(output_digits<REAL>());
    std::cout << "  max radius = " << max_radius << '\n';
}

} // namespace

int main() {
    show_tier<float>("float");
    show_tier<double>("double");
    mpfrxx::set_default_precision_bits(512);
    show_tier<mpfrxx::mpfr_class>("mpfr@512");
    return EXIT_SUCCESS;
}
