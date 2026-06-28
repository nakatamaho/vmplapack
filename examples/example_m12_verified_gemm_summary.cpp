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

constexpr std::ptrdiff_t rows = 4;
constexpr std::ptrdiff_t inner = 10;
constexpr std::ptrdiff_t cols = 4;

template <class REAL>
int output_digits() {
    return std::numeric_limits<REAL>::max_digits10;
}

template <>
int output_digits<mpfrxx::mpfr_class>() {
    return 24;
}

int oracle_output_digits() {
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
REAL ratio_value(int numerator, int denominator) {
    REAL n = from_int<REAL>(numerator);
    REAL d = from_int<REAL>(denominator);
    REAL q = n / d;
    return q;
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
mpfrxx::mpfr_class midpoint_error(const vmplapack::Rmidrad<REAL>& box,
                                  const mpfrxx::mpfr_class& oracle_mid,
                                  mpfr_prec_t precision) {
    mpfrxx::mpfr_class mid = vmplapack::oracle::widen_value(box.mid, precision);
    mpfrxx::mpfr_class ref = vmplapack::oracle::widen_mpfr(oracle_mid, precision);
    mpfrxx::mpfr_class diff = mpfrxx::mpfr_class::with_precision(precision);
    diff = mid - ref;
    mpfrxx::mpfr_class abs_diff = vmplapack::oracle::abs_at(diff, precision);
    return abs_diff;
}

template <class REAL>
bool midpoint_differs(const vmplapack::Rmidrad<REAL>& box,
                      const mpfrxx::mpfr_class& oracle_mid,
                      mpfr_prec_t precision) {
    mpfrxx::mpfr_class mid = vmplapack::oracle::widen_value(box.mid, precision);
    mpfrxx::mpfr_class ref = vmplapack::oracle::widen_mpfr(oracle_mid, precision);
    return mpfr_cmp(mid.mpfr_data(), ref.mpfr_data()) != 0;
}

template <class REAL>
void append_problem(std::vector<REAL>& left, std::vector<REAL>& right, int exponent) {
    using A = vmplapack::Rarith<REAL>;

    REAL big = power_of_two_value<REAL>(exponent);
    REAL one = A::one();
    REAL minus_one = -one;

    for (std::ptrdiff_t row = 0; row < rows; ++row) {
        left.push_back(big);
        left.push_back(-big);
        for (std::ptrdiff_t k = 2; k < inner; ++k) {
            int numerator = static_cast<int>((row + 2) * (k + 3) + 1);
            int denominator = static_cast<int>((row + 3) * (k + 5) + 7);
            REAL value = ratio_value<REAL>(numerator, denominator);
            if (((row + k) % 2) != 0) {
                value = -value;
            }
            left.push_back(value);
        }
    }

    for (std::ptrdiff_t k = 0; k < inner; ++k) {
        for (std::ptrdiff_t col = 0; col < cols; ++col) {
            if (k < 2) {
                right.push_back(one);
            } else {
                int numerator = static_cast<int>((col + 3) * (k + 2) + 2);
                int denominator = static_cast<int>((col + 5) * (k + 4) + 11);
                REAL value = ratio_value<REAL>(numerator, denominator);
                if (((col + 2 * k) % 3) == 0) {
                    value = minus_one * value;
                }
                right.push_back(value);
            }
        }
    }
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
            std::cout << "]; ";
        } else {
            std::cout << "] ";
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
            std::cout << "]; ";
        } else {
            std::cout << "] ";
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
            std::cout << "]; ";
        } else {
            std::cout << "] ";
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
            std::cout << "]; ";
        } else {
            std::cout << "] ";
        }
    }
    std::cout << "]";
}

template <class REAL>
void compute_oracle_matrix(const std::vector<REAL>& left,
                           const std::vector<REAL>& right,
                           std::vector<vmplapack::oracle::Rdot_interval>& intervals,
                           std::vector<mpfrxx::mpfr_class>& oracle_midpoints) {
    intervals.clear();
    oracle_midpoints.clear();
    intervals.reserve(static_cast<std::size_t>(rows * cols));
    oracle_midpoints.reserve(static_cast<std::size_t>(rows * cols));
    for (std::ptrdiff_t row = 0; row < rows; ++row) {
        for (std::ptrdiff_t col = 0; col < cols; ++col) {
            vmplapack::oracle::Rdot_interval ref =
                vmplapack::oracle::Rdot_oracle(inner, left.data() + row * inner, 1, right.data() + col, cols);
            mpfr_prec_t precision = ref.precision + static_cast<mpfr_prec_t>(vmplapack::Rarith<REAL>::precision_bits() + 128L);
            mpfrxx::mpfr_class mid = interval_midpoint(ref, precision);
            intervals.push_back(ref);
            oracle_midpoints.push_back(mid);
        }
    }
}

template <class REAL>
void show_tier(const char* tier, int exponent) {
    using A = vmplapack::Rarith<REAL>;

    std::vector<REAL> left;
    std::vector<REAL> right;
    append_problem(left, right, exponent);

    std::vector<vmplapack::Rmidrad<REAL>> result(static_cast<std::size_t>(rows * cols));
    vmplapack::VerificationStatus status =
        vmplapack::vRgemm_point<REAL>(rows, cols, inner, left.data(), inner, right.data(), cols, result.data(), cols);

    std::vector<vmplapack::oracle::Rdot_interval> intervals;
    std::vector<mpfrxx::mpfr_class> oracle_midpoints;
    compute_oracle_matrix(left, right, intervals, oracle_midpoints);

    bool all_ok = true;
    bool all_covered = true;
    int midpoint_mismatch_count = 0;
    REAL max_radius = A::zero();
    mpfr_prec_t error_precision = static_cast<mpfr_prec_t>(std::max(1024L, 4L * A::precision_bits() + 256L));
    mpfrxx::mpfr_class max_midpoint_error = vmplapack::oracle::zero_at(error_precision);
    std::ptrdiff_t worst_index = 0;
    for (std::ptrdiff_t i = 0; i < rows * cols; ++i) {
        const vmplapack::Rmidrad<REAL>& box = result[static_cast<std::size_t>(i)];
        const vmplapack::oracle::Rdot_interval& ref = intervals[static_cast<std::size_t>(i)];
        const mpfrxx::mpfr_class& oracle_mid = oracle_midpoints[static_cast<std::size_t>(i)];
        if (box.status != vmplapack::Rstatus::ok) {
            all_ok = false;
        }
        if (!midrad_covers_oracle(box, ref)) {
            all_covered = false;
        }
        if (midpoint_differs(box, oracle_mid, error_precision)) {
            ++midpoint_mismatch_count;
        }
        mpfrxx::mpfr_class error = midpoint_error(box, oracle_mid, error_precision);
        if (error > max_midpoint_error) {
            max_midpoint_error = error;
        }
        if (box.rad > max_radius) {
            max_radius = box.rad;
            worst_index = i;
        }
    }

    std::ptrdiff_t worst_row = worst_index / cols;
    std::ptrdiff_t worst_col = worst_index % cols;
    const vmplapack::Rmidrad<REAL>& worst = result[static_cast<std::size_t>(worst_index)];
    REAL worst_lo = lower_display(worst);
    REAL worst_hi = upper_display(worst);

    std::cout << std::setprecision(output_digits<REAL>()) << std::boolalpha;
    std::cout << tier << '\n';
    std::cout << "  problem = A(4x10) * B(10x4), row-major" << '\n';
    std::cout << "  cancellation scale = 2^" << exponent << '\n';
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
    std::cout << "  return status = " << status_name(status) << '\n';
    std::cout << "  all components ok = " << all_ok << '\n';
    std::cout << "  all oracle intervals covered = " << all_covered << '\n';
    std::cout << "  components with C.mid != oracle midpoint = " << midpoint_mismatch_count << " / " << rows * cols << '\n';
    std::cout << std::setprecision(oracle_output_digits());
    std::cout << "  max |C.mid - oracle midpoint| = " << max_midpoint_error << '\n';
    std::cout << std::setprecision(output_digits<REAL>());
    std::cout << "  max radius = " << max_radius << '\n';
    std::cout << "  worst radius component = C(" << worst_row << "," << worst_col << ")" << '\n';
    std::cout << "    oracle midpoint = " << std::setprecision(oracle_output_digits())
              << oracle_midpoints[static_cast<std::size_t>(worst_index)] << '\n';
    std::cout << std::setprecision(output_digits<REAL>());
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
