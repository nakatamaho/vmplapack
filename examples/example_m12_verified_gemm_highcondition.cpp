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
    int float_cond_power;
    int double_cond_power;
    int mpfr_cond_power;
    std::uint64_t seed;
    bool no_matrices;
};

struct tier_summary {
    bool ok;
    const char* message;
};

constexpr unsigned dyadic_random_bits = 10U;
constexpr unsigned dyadic_random_bins = 1U << dyadic_random_bits;
constexpr unsigned residual_numerator_min = 65U;
constexpr unsigned residual_numerator_span = 63U;
constexpr unsigned residual_denominator_min = 129U;
constexpr unsigned residual_denominator_span = 127U;

options default_options() {
    options opt;
    opt.m = 4;
    opt.n = 4;
    opt.k = 9;
    opt.float_cond_power = 24;
    opt.double_cond_power = 52;
    opt.mpfr_cond_power = 256;
    opt.seed = 17U;
    opt.no_matrices = false;
    return opt;
}

void print_usage(const char* program) {
    std::cout << "usage: " << program << " [options]\n"
              << "  --m M                  output rows (default 4)\n"
              << "  --n N                  output columns (default 4)\n"
              << "  --k K                  inner dimension, K >= 3 (default 9)\n"
              << "  --cond P               set all target ORO condition powers to 2^P\n"
              << "  --float-cond P         float target ORO condition power (default 24)\n"
              << "  --double-cond P        double target ORO condition power (default 52)\n"
              << "  --mpfr-cond P          mpfr@512 target ORO condition power (default 256)\n"
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
        } else if (std::strcmp(argv[i], "--cond") == 0) {
            if (!parse_long(argv[++i], &parsed)) {
                std::cerr << "invalid --cond\n";
                return false;
            }
            opt->float_cond_power = static_cast<int>(parsed);
            opt->double_cond_power = static_cast<int>(parsed);
            opt->mpfr_cond_power = static_cast<int>(parsed);
        } else if (std::strcmp(argv[i], "--float-cond") == 0) {
            if (!parse_long(argv[++i], &parsed)) {
                std::cerr << "invalid --float-cond\n";
                return false;
            }
            opt->float_cond_power = static_cast<int>(parsed);
        } else if (std::strcmp(argv[i], "--double-cond") == 0) {
            if (!parse_long(argv[++i], &parsed)) {
                std::cerr << "invalid --double-cond\n";
                return false;
            }
            opt->double_cond_power = static_cast<int>(parsed);
        } else if (std::strcmp(argv[i], "--mpfr-cond") == 0) {
            if (!parse_long(argv[++i], &parsed)) {
                std::cerr << "invalid --mpfr-cond\n";
                return false;
            }
            opt->mpfr_cond_power = static_cast<int>(parsed);
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
    if (opt.m <= 0 || opt.n <= 0 || opt.k < 3) {
        std::cerr << "requires --m > 0, --n > 0, and --k >= 3\n";
        return false;
    }
    if (opt.float_cond_power < 1 || opt.double_cond_power < 1 || opt.mpfr_cond_power < 1) {
        std::cerr << "condition powers must be positive\n";
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
        static_assert(vmplapack::Ralways_false<REAL>::value, "max_power_exponent does not support this REAL type.");
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
        static_assert(vmplapack::Ralways_false<REAL>::value, "min_power_exponent does not support this REAL type.");
    }
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
REAL random_dyadic_unit(std::mt19937_64* engine) {
    unsigned bucket = static_cast<unsigned>((*engine)() & static_cast<std::uint64_t>(dyadic_random_bins - 1U));
    REAL numerator = from_int<REAL>(static_cast<int>(dyadic_random_bins + bucket));
    REAL denominator = from_int<REAL>(static_cast<int>(2U * dyadic_random_bins));
    REAL value = numerator / denominator;
    return value;
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
REAL random_signed_dyadic_unit(std::mt19937_64* engine) {
    REAL scale = random_dyadic_unit<REAL>(engine);
    REAL sign = random_sign<REAL>(engine);
    REAL value = sign * scale;
    return value;
}

template <class REAL>
REAL random_signed_residual(std::mt19937_64* engine) {
    unsigned raw_numerator =
        static_cast<unsigned>((*engine)() % static_cast<std::uint64_t>(residual_numerator_span)) + residual_numerator_min;
    unsigned raw_denominator =
        static_cast<unsigned>((*engine)() % static_cast<std::uint64_t>(residual_denominator_span)) + residual_denominator_min;
    int numerator = static_cast<int>(raw_numerator);
    int denominator = static_cast<int>(raw_denominator | 1U);
    REAL magnitude = ratio_value<REAL>(numerator, denominator);
    REAL sign = random_sign<REAL>(engine);
    REAL value = sign * magnitude;
    return value;
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
tier_summary build_high_condition_gemm(std::ptrdiff_t m,
                                       std::ptrdiff_t n,
                                       std::ptrdiff_t k,
                                       int cond_power,
                                       std::uint64_t seed,
                                       std::vector<REAL>& left,
                                       std::vector<REAL>& right,
                                       int* scale_out,
                                       std::ptrdiff_t* pairs_out) {
    using A = vmplapack::Rarith<REAL>;

    std::ptrdiff_t pairs = (k - 1) / 2;
    if (pairs < 1) {
        return {false, "k is too small for a high-condition construction"};
    }

    int scale = cond_power - ceil_log2_size(static_cast<std::size_t>(2 * pairs));
    if (scale < min_power_exponent<REAL>() || scale > max_power_exponent<REAL>()) {
        return {false, "requested condition is not representable for this tier"};
    }

    left.assign(static_cast<std::size_t>(m * k), A::zero());
    right.assign(static_cast<std::size_t>(k * n), A::zero());

    std::mt19937_64 engine(seed);
    REAL big = power_of_two_value<REAL>(scale);
    std::ptrdiff_t residual_index = 2 * pairs;

    std::vector<REAL> row_residual(static_cast<std::size_t>(m));
    std::vector<REAL> col_residual(static_cast<std::size_t>(n));
    for (std::ptrdiff_t row = 0; row < m; ++row) {
        row_residual[static_cast<std::size_t>(row)] = random_signed_residual<REAL>(&engine);
    }
    for (std::ptrdiff_t col = 0; col < n; ++col) {
        col_residual[static_cast<std::size_t>(col)] = random_signed_residual<REAL>(&engine);
    }

    for (std::ptrdiff_t row = 0; row < m; ++row) {
        for (std::ptrdiff_t pair = 0; pair < pairs; ++pair) {
            REAL pair_scale = random_dyadic_unit<REAL>(&engine);
            REAL high = row_residual[static_cast<std::size_t>(row)] * big;
            REAL high_scaled = high * pair_scale;
            left[static_cast<std::size_t>(row * k + 2 * pair)] = high_scaled;
            left[static_cast<std::size_t>(row * k + 2 * pair + 1)] = high_scaled;
        }
        left[static_cast<std::size_t>(row * k + residual_index)] = row_residual[static_cast<std::size_t>(row)];
    }

    for (std::ptrdiff_t col = 0; col < n; ++col) {
        for (std::ptrdiff_t pair = 0; pair < pairs; ++pair) {
            REAL pair_scale = random_dyadic_unit<REAL>(&engine);
            REAL high_col = col_residual[static_cast<std::size_t>(col)] * pair_scale;
            REAL sign = random_sign<REAL>(&engine);
            REAL first = high_col * sign;
            REAL second = -first;
            right[static_cast<std::size_t>((2 * pair) * n + col)] = first;
            right[static_cast<std::size_t>((2 * pair + 1) * n + col)] = second;
        }
        right[static_cast<std::size_t>(residual_index * n + col)] = col_residual[static_cast<std::size_t>(col)];
    }

    *scale_out = scale;
    *pairs_out = pairs;
    return {true, "ok"};
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

bool should_print_matrices(const options& opt) {
    return !opt.no_matrices;
}

template <class REAL>
void show_tier(const char* tier,
               std::ptrdiff_t m,
               std::ptrdiff_t n,
               std::ptrdiff_t k,
               int cond_power,
               std::uint64_t seed,
               bool print_matrices) {
    using A = vmplapack::Rarith<REAL>;

    std::vector<REAL> left;
    std::vector<REAL> right;
    int scale = 0;
    std::ptrdiff_t pairs = 0;
    tier_summary generated = build_high_condition_gemm(m, n, k, cond_power, seed, left, right, &scale, &pairs);

    std::cout << std::setprecision(output_digits<REAL>()) << std::boolalpha;
    std::cout << tier << '\n';
    std::cout << "  problem = A(" << m << "x" << k << ") * B(" << k << "x" << n << "), row-major" << '\n';
    std::cout << "  target cond_oro = 2^" << cond_power << '\n';
    std::cout << "  seed = " << static_cast<unsigned long long>(seed) << '\n';
    std::cout << "  generator status = " << generated.message << '\n';
    if (!generated.ok) {
        return;
    }
    std::cout << "  construction = randomized dyadic cancellation pairs plus one residual term" << '\n';
    std::cout << "  cancellation pairs per component = " << pairs << '\n';
    std::cout << "  cancellation scale = 2^" << scale << '\n';
    std::cout << "  dyadic scale bits = " << dyadic_random_bits << '\n';
    std::cout << "  residual numerator range = [" << residual_numerator_min << ", "
              << (residual_numerator_min + residual_numerator_span - 1U) << "]" << '\n';
    std::cout << "  residual odd denominator range = [" << residual_denominator_min << ", "
              << (residual_denominator_min + residual_denominator_span - 1U) << "]" << '\n';

    std::vector<vmplapack::Rmidrad<REAL>> result(static_cast<std::size_t>(m * n));
    vmplapack::VerificationStatus status =
        vmplapack::vRgemm_point<REAL>(m, n, k, left.data(), k, right.data(), n, result.data(), n);

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
    oracle_midpoints.reserve(static_cast<std::size_t>(m * n));
    std::vector<mpfrxx::mpfr_class> midpoint_diffs;
    midpoint_diffs.reserve(static_cast<std::size_t>(m * n));

    for (std::ptrdiff_t row = 0; row < m; ++row) {
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            std::ptrdiff_t index = row * n + col;
            const vmplapack::Rmidrad<REAL>& box = result[static_cast<std::size_t>(index)];
            vmplapack::oracle::Rdot_interval ref =
                vmplapack::oracle::Rdot_oracle(k, left.data() + row * k, 1, right.data() + col, n);
            mpfr_prec_t local_precision = ref.precision + static_cast<mpfr_prec_t>(A::precision_bits() + 128L);
            mpfrxx::mpfr_class oracle_mid = interval_midpoint(ref, local_precision);
            oracle_midpoints.push_back(oracle_mid);
            mpfrxx::mpfr_class sum_abs =
                vmplapack::oracle::upward_abs_term_sum(k, left.data() + row * k, 1, right.data() + col, n, local_precision);
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

    std::ptrdiff_t worst_row = worst_index / n;
    std::ptrdiff_t worst_col = worst_index % n;
    const vmplapack::Rmidrad<REAL>& worst = result[static_cast<std::size_t>(worst_index)];
    REAL worst_lo = lower_display(worst);
    REAL worst_hi = upper_display(worst);

    if (print_matrices) {
        std::cout << "  A = ";
        print_octave_matrix(m, k, left.data(), k);
        std::cout << '\n';
        std::cout << "  B = ";
        print_octave_matrix(k, n, right.data(), n);
        std::cout << '\n';
        std::cout << "  oracle C midpoint = ";
        std::cout << std::setprecision(oracle_output_digits());
        print_mpfr_matrix(m, n, oracle_midpoints, n);
        std::cout << '\n';
        std::cout << std::setprecision(output_digits<REAL>());
        std::cout << "  result C.mid = ";
        print_mid_matrix(m, n, result, n);
        std::cout << '\n';
        std::cout << "  result C.rad = ";
        print_rad_matrix(m, n, result, n);
        std::cout << '\n';
        std::cout << "  result C.diff = ";
        std::cout << std::scientific << std::setprecision(diff_output_digits());
        print_mpfr_matrix(m, n, midpoint_diffs, n);
        std::cout << '\n';
        std::cout << std::defaultfloat << std::setprecision(output_digits<REAL>());
    }

    std::cout << "  return status = " << status_name(status) << '\n';
    std::cout << "  all components ok = " << all_ok << '\n';
    std::cout << "  all oracle intervals covered = " << all_covered << '\n';
    std::cout << std::setprecision(oracle_output_digits());
    std::cout << "  min measured cond_oro = " << min_condition << '\n';
    std::cout << "  max measured cond_oro = " << max_condition << '\n';
    std::cout << "  components with C.mid != oracle midpoint = " << midpoint_mismatch_count << " / " << m * n << '\n';
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

    bool print_matrices = should_print_matrices(opt);
    show_tier<float>("float", opt.m, opt.n, opt.k, opt.float_cond_power, opt.seed, print_matrices);
    show_tier<double>("double", opt.m, opt.n, opt.k, opt.double_cond_power, opt.seed + 1009U, print_matrices);
    mpfrxx::set_default_precision_bits(512);
    show_tier<mpfrxx::mpfr_class>("mpfr@512", opt.m, opt.n, opt.k, opt.mpfr_cond_power, opt.seed + 2003U, print_matrices);
    return EXIT_SUCCESS;
}
