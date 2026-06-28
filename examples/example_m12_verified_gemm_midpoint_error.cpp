// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include "../tests/Rdot_oracle.h"

#include <vmplapack/vmplapack.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <type_traits>
#include <vector>

namespace {

struct options {
    std::ptrdiff_t m;
    std::ptrdiff_t n;
    std::ptrdiff_t k;
    std::uint64_t seed;
    bool no_matrices;
};

constexpr unsigned rational_numerator_min = 1U;
constexpr unsigned rational_numerator_span = 127U;
constexpr unsigned rational_denominator_min = 65U;
constexpr unsigned rational_denominator_span = 191U;

options default_options() {
    options opt;
    opt.m = 4;
    opt.n = 4;
    opt.k = 10;
    opt.seed = 17U;
    opt.no_matrices = false;
    return opt;
}

void print_usage(const char* program) {
    std::cout << "usage: " << program << " [options]\n"
              << "  --m M                  output rows (default 4)\n"
              << "  --n N                  output columns (default 4)\n"
              << "  --k K                  inner dimension, K > 0 (default 10)\n"
              << "  --seed S               random seed (default 17)\n"
              << "  --no-matrices          suppress A, B, C.mid, C.rad, C.diff\n";
}

bool parse_long(const char* text, long* value) {
    char* end = nullptr;
    long parsed = std::strtol(text, &end, 10);
    if (end == text || *end != '\0') {
        return false;
    }
    *value = parsed;
    return true;
}

bool parse_u64(const char* text, std::uint64_t* value) {
    char* end = nullptr;
    unsigned long long parsed = std::strtoull(text, &end, 10);
    if (end == text || *end != '\0') {
        return false;
    }
    *value = static_cast<std::uint64_t>(parsed);
    return true;
}

bool parse_options(int argc, char** argv, options* opt) {
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            std::exit(EXIT_SUCCESS);
        }
        if (std::strcmp(argv[i], "--no-matrices") == 0) {
            opt->no_matrices = true;
            continue;
        }

        if (i + 1 >= argc) {
            std::cerr << "missing value after " << argv[i] << '\n';
            return false;
        }

        long parsed = 0;
        if (std::strcmp(argv[i], "--m") == 0) {
            if (!parse_long(argv[++i], &parsed)) {
                std::cerr << "invalid --m\n";
                return false;
            }
            opt->m = static_cast<std::ptrdiff_t>(parsed);
        } else if (std::strcmp(argv[i], "--n") == 0) {
            if (!parse_long(argv[++i], &parsed)) {
                std::cerr << "invalid --n\n";
                return false;
            }
            opt->n = static_cast<std::ptrdiff_t>(parsed);
        } else if (std::strcmp(argv[i], "--k") == 0) {
            if (!parse_long(argv[++i], &parsed)) {
                std::cerr << "invalid --k\n";
                return false;
            }
            opt->k = static_cast<std::ptrdiff_t>(parsed);
        } else if (std::strcmp(argv[i], "--seed") == 0) {
            std::uint64_t parsed_seed = 0U;
            if (!parse_u64(argv[++i], &parsed_seed)) {
                std::cerr << "invalid --seed\n";
                return false;
            }
            opt->seed = parsed_seed;
        } else {
            std::cerr << "unknown option: " << argv[i] << '\n';
            return false;
        }
    }
    return true;
}

bool validate_options(const options& opt) {
    if (opt.m <= 0 || opt.n <= 0 || opt.k <= 0) {
        std::cerr << "requires --m > 0, --n > 0, and --k > 0\n";
        return false;
    }
    return true;
}

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
REAL ratio_value(int numerator, int denominator) {
    REAL n = from_int<REAL>(numerator);
    REAL d = from_int<REAL>(denominator);
    REAL q = n / d;
    return q;
}

template <class REAL>
REAL random_sign(std::mt19937_64* engine) {
    REAL one = vmplapack::Rarith<REAL>::one();
    if (((*engine)() & 1ULL) == 0ULL) {
        return one;
    }
    REAL minus_one = -one;
    return minus_one;
}

template <class REAL>
REAL random_rational_value(std::mt19937_64* engine) {
    unsigned raw_numerator =
        static_cast<unsigned>((*engine)() % static_cast<std::uint64_t>(rational_numerator_span)) + rational_numerator_min;
    unsigned raw_denominator =
        static_cast<unsigned>((*engine)() % static_cast<std::uint64_t>(rational_denominator_span)) + rational_denominator_min;
    int numerator = static_cast<int>(raw_numerator);
    int denominator = static_cast<int>(raw_denominator | 1U);
    REAL value = ratio_value<REAL>(numerator, denominator);
    REAL sign = random_sign<REAL>(engine);
    REAL signed_value = sign * value;
    return signed_value;
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
bool midpoint_differs(const vmplapack::Rmidrad<REAL>& box,
                      const mpfrxx::mpfr_class& oracle_mid,
                      mpfr_prec_t precision) {
    mpfrxx::mpfr_class mid = vmplapack::oracle::widen_value(box.mid, precision);
    mpfrxx::mpfr_class ref = vmplapack::oracle::widen_mpfr(oracle_mid, precision);
    return mpfr_cmp(mid.mpfr_data(), ref.mpfr_data()) != 0;
}

mpfrxx::mpfr_class condition_from_oracle_mid(const mpfrxx::mpfr_class& oracle_mid,
                                             const mpfrxx::mpfr_class& sum_abs,
                                             mpfr_prec_t precision) {
    mpfrxx::mpfr_class abs_mid = vmplapack::oracle::abs_at(oracle_mid, precision);
    if (abs_mid == vmplapack::oracle::zero_at(precision)) {
        mpfrxx::mpfr_class inf = mpfrxx::mpfr_class::with_precision(precision);
        mpfr_set_inf(inf.mpfr_data(), 1);
        return inf;
    }
    mpfrxx::mpfr_class two = mpfrxx::mpfr_class::with_precision(precision, 2.0);
    mpfrxx::mpfr_class numerator = mpfrxx::mpfr_class::with_precision(precision);
    numerator = two * sum_abs;
    mpfrxx::mpfr_class condition = mpfrxx::mpfr_class::with_precision(precision);
    {
        mpfrxx::rounding_mode_scope scope(MPFR_RNDU);
        condition = numerator / abs_mid;
    }
    return condition;
}

template <class REAL>
void build_random_midpoint_problem(std::ptrdiff_t m,
                                   std::ptrdiff_t n,
                                   std::ptrdiff_t k,
                                   std::uint64_t seed,
                                   std::vector<REAL>& left,
                                   std::vector<REAL>& right) {
    std::mt19937_64 engine(seed);
    left.clear();
    right.clear();
    left.reserve(static_cast<std::size_t>(m * k));
    right.reserve(static_cast<std::size_t>(k * n));

    for (std::ptrdiff_t row = 0; row < m; ++row) {
        for (std::ptrdiff_t pos = 0; pos < k; ++pos) {
            REAL value = random_rational_value<REAL>(&engine);
            left.push_back(value);
        }
    }
    for (std::ptrdiff_t pos = 0; pos < k; ++pos) {
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            REAL value = random_rational_value<REAL>(&engine);
            right.push_back(value);
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
void show_tier(const char* tier,
               const options& opt,
               std::uint64_t tier_seed) {
    using A = vmplapack::Rarith<REAL>;

    std::vector<REAL> left;
    std::vector<REAL> right;
    build_random_midpoint_problem<REAL>(opt.m, opt.n, opt.k, tier_seed, left, right);

    std::vector<vmplapack::Rmidrad<REAL>> result(static_cast<std::size_t>(opt.m * opt.n));
    vmplapack::VerificationStatus status =
        vmplapack::vRgemm_point<REAL>(opt.m, opt.n, opt.k, left.data(), opt.k, right.data(), opt.n, result.data(), opt.n);

    bool all_ok = true;
    bool all_covered = true;
    int midpoint_mismatch_count = 0;
    REAL max_radius = A::zero();
    std::ptrdiff_t worst_index = 0;
    mpfr_prec_t metric_precision = static_cast<mpfr_prec_t>(std::max(1024L, 4L * A::precision_bits() + 256L));
    mpfrxx::mpfr_class max_midpoint_error = vmplapack::oracle::zero_at(metric_precision);
    mpfrxx::mpfr_class min_condition = mpfrxx::mpfr_class::with_precision(metric_precision);
    mpfrxx::mpfr_class max_condition = vmplapack::oracle::zero_at(metric_precision);
    bool have_condition = false;
    std::vector<mpfrxx::mpfr_class> oracle_midpoints;
    std::vector<mpfrxx::mpfr_class> midpoint_diffs;
    oracle_midpoints.reserve(static_cast<std::size_t>(opt.m * opt.n));
    midpoint_diffs.reserve(static_cast<std::size_t>(opt.m * opt.n));

    for (std::ptrdiff_t row = 0; row < opt.m; ++row) {
        for (std::ptrdiff_t col = 0; col < opt.n; ++col) {
            std::ptrdiff_t index = row * opt.n + col;
            const vmplapack::Rmidrad<REAL>& box = result[static_cast<std::size_t>(index)];
            vmplapack::oracle::Rdot_interval ref =
                vmplapack::oracle::Rdot_oracle(opt.k, left.data() + row * opt.k, 1, right.data() + col, opt.n);
            mpfr_prec_t local_precision = ref.precision + static_cast<mpfr_prec_t>(A::precision_bits() + 128L);
            mpfrxx::mpfr_class oracle_mid = interval_midpoint(ref, local_precision);
            oracle_midpoints.push_back(oracle_mid);
            mpfrxx::mpfr_class sum_abs =
                vmplapack::oracle::upward_abs_term_sum(opt.k, left.data() + row * opt.k, 1, right.data() + col, opt.n,
                                                       local_precision);
            mpfrxx::mpfr_class condition = condition_from_oracle_mid(oracle_mid, sum_abs, local_precision);
            mpfrxx::mpfr_class condition_metric = vmplapack::oracle::widen_mpfr(condition, metric_precision);
            if (!have_condition || condition_metric < min_condition) {
                min_condition = condition_metric;
            }
            if (condition_metric > max_condition) {
                max_condition = condition_metric;
            }
            have_condition = true;

            if (box.status != vmplapack::Rstatus::ok) {
                all_ok = false;
            }
            if (!midrad_covers_oracle(box, ref)) {
                all_covered = false;
            }
            if (midpoint_differs(box, oracle_mid, metric_precision)) {
                ++midpoint_mismatch_count;
            }
            midpoint_diffs.push_back(midpoint_difference(box, oracle_mid, metric_precision));
            mpfrxx::mpfr_class error = midpoint_error(box, oracle_mid, metric_precision);
            if (error > max_midpoint_error) {
                max_midpoint_error = error;
            }
            if (box.rad > max_radius) {
                max_radius = box.rad;
                worst_index = index;
            }
        }
    }

    std::ptrdiff_t worst_row = worst_index / opt.n;
    std::ptrdiff_t worst_col = worst_index % opt.n;
    const vmplapack::Rmidrad<REAL>& worst = result[static_cast<std::size_t>(worst_index)];
    REAL worst_lo = lower_display(worst);
    REAL worst_hi = upper_display(worst);

    std::cout << std::setprecision(output_digits<REAL>()) << std::boolalpha;
    std::cout << tier << '\n';
    std::cout << "  problem = A(" << opt.m << "x" << opt.k << ") * B(" << opt.k << "x" << opt.n << "), row-major" << '\n';
    std::cout << "  construction = independent seeded signed rational entries" << '\n';
    std::cout << "  seed = " << static_cast<unsigned long long>(tier_seed) << '\n';
    std::cout << "  numerator range = [" << rational_numerator_min << ", "
              << (rational_numerator_min + rational_numerator_span - 1U) << "]" << '\n';
    std::cout << "  odd denominator range = [" << rational_denominator_min << ", "
              << (rational_denominator_min + rational_denominator_span - 1U) << "]" << '\n';
    if (!opt.no_matrices) {
        std::cout << "  A = ";
        print_octave_matrix(opt.m, opt.k, left.data(), opt.k);
        std::cout << '\n';
        std::cout << "  B = ";
        print_octave_matrix(opt.k, opt.n, right.data(), opt.n);
        std::cout << '\n';
        std::cout << "  oracle C midpoint = ";
        std::cout << std::setprecision(oracle_output_digits());
        print_mpfr_matrix(opt.m, opt.n, oracle_midpoints, opt.n);
        std::cout << '\n';
        std::cout << std::setprecision(output_digits<REAL>());
        std::cout << "  result C.mid = ";
        print_mid_matrix(opt.m, opt.n, result, opt.n);
        std::cout << '\n';
        std::cout << "  result C.rad = ";
        print_rad_matrix(opt.m, opt.n, result, opt.n);
        std::cout << '\n';
        std::cout << "  result C.diff = ";
        std::cout << std::scientific << std::setprecision(diff_output_digits());
        print_mpfr_matrix(opt.m, opt.n, midpoint_diffs, opt.n);
        std::cout << '\n';
        std::cout << std::defaultfloat << std::setprecision(output_digits<REAL>());
    }
    std::cout << "  return status = " << status_name(status) << '\n';
    std::cout << "  all components ok = " << all_ok << '\n';
    std::cout << "  all oracle intervals covered = " << all_covered << '\n';
    std::cout << std::setprecision(oracle_output_digits());
    std::cout << "  min measured cond_oro = " << min_condition << '\n';
    std::cout << "  max measured cond_oro = " << max_condition << '\n';
    std::cout << "  components with C.mid != oracle midpoint = " << midpoint_mismatch_count << " / " << opt.m * opt.n << '\n';
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

int main(int argc, char** argv) {
    options opt = default_options();
    if (!parse_options(argc, argv, &opt) || !validate_options(opt)) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    show_tier<float>("float", opt, opt.seed);
    show_tier<double>("double", opt, opt.seed + 1009U);
    mpfrxx::set_default_precision_bits(512);
    show_tier<mpfrxx::mpfr_class>("mpfr@512", opt, opt.seed + 2003U);
    return EXIT_SUCCESS;
}
