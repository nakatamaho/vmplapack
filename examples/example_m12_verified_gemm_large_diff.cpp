// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include "../tests/Rdot_oracle.h"

#include <vmplapack/vmplapack.h>

#include <algorithm>
#include <cmath>
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

constexpr unsigned dyadic_random_bits = 10U;
constexpr unsigned dyadic_random_bins = 1U << dyadic_random_bits;

options default_options() {
    options opt;
    opt.m = 3;
    opt.n = 3;
    opt.k = 9;
    opt.seed = 29U;
    opt.no_matrices = false;
    return opt;
}

void print_usage(const char* program) {
    std::cout << "usage: " << program << " [options]\n"
              << "  --m M                  output rows (default 3)\n"
              << "  --n N                  output columns (default 3)\n"
              << "  --k K                  inner dimension, K >= 4 (default 9)\n"
              << "  --seed S               random seed (default 29)\n"
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
    if (opt.m <= 0 || opt.n <= 0 || opt.k < 4) {
        std::cerr << "requires --m > 0, --n > 0, and --k >= 4\n";
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
REAL from_unsigned(unsigned value) {
    REAL result = static_cast<REAL>(value);
    return result;
}

template <>
mpfrxx::mpfr_class from_unsigned<mpfrxx::mpfr_class>(unsigned value) {
    mpfrxx::mpfr_class result =
        mpfrxx::mpfr_class::with_precision(static_cast<mpfr_prec_t>(vmplapack::Rarith<mpfrxx::mpfr_class>::precision_bits()));
    mpfr_set_ui(result.mpfr_data(), static_cast<unsigned long>(value), MPFR_RNDN);
    return result;
}

template <class REAL>
REAL random_dyadic_unit(std::mt19937_64* engine) {
    unsigned bucket = static_cast<unsigned>((*engine)() & static_cast<std::uint64_t>(dyadic_random_bins - 1U));
    REAL numerator = from_unsigned<REAL>(dyadic_random_bins + bucket);
    REAL denominator = from_unsigned<REAL>(2U * dyadic_random_bins);
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
REAL signed_scaled_power(int exponent, REAL sign, REAL scale) {
    REAL base = power_of_two_value<REAL>(exponent);
    REAL scaled = base * scale;
    REAL value = sign * scaled;
    return value;
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
void build_random_large_diff_problem(std::ptrdiff_t m,
                                     std::ptrdiff_t n,
                                     std::ptrdiff_t k,
                                     std::uint64_t seed,
                                     std::vector<REAL>& left,
                                     std::vector<REAL>& right) {
    using A = vmplapack::Rarith<REAL>;

    std::mt19937_64 engine(seed);
    left.assign(static_cast<std::size_t>(m * k), A::zero());
    right.assign(static_cast<std::size_t>(k * n), A::zero());

    std::vector<REAL> h1_row_sign(static_cast<std::size_t>(m));
    std::vector<REAL> h2_row_sign(static_cast<std::size_t>(m));
    std::vector<REAL> h1_col_sign(static_cast<std::size_t>(n));
    std::vector<REAL> h2_col_sign(static_cast<std::size_t>(n));
    for (std::ptrdiff_t row = 0; row < m; ++row) {
        h1_row_sign[static_cast<std::size_t>(row)] = random_sign<REAL>(&engine);
        h2_row_sign[static_cast<std::size_t>(row)] = random_sign<REAL>(&engine);
    }
    for (std::ptrdiff_t col = 0; col < n; ++col) {
        h1_col_sign[static_cast<std::size_t>(col)] = random_sign<REAL>(&engine);
        h2_col_sign[static_cast<std::size_t>(col)] = random_sign<REAL>(&engine);
    }

    std::ptrdiff_t blocks = k / 4;
    for (std::ptrdiff_t block = 0; block < blocks; ++block) {
        std::ptrdiff_t base = 4 * block;
        for (std::ptrdiff_t row = 0; row < m; ++row) {
            REAL p_scale = random_dyadic_unit<REAL>(&engine);
            REAL p_sign = random_sign<REAL>(&engine);
            REAL p_value = signed_scaled_power<REAL>(large_diff_exponents<REAL>::p, p_sign, p_scale);
            REAL neg_p_value = -p_value;
            left[static_cast<std::size_t>(row * k + base)] = p_value;
            left[static_cast<std::size_t>(row * k + base + 3)] = neg_p_value;

            REAL h1_scale = random_dyadic_unit<REAL>(&engine);
            REAL h1_sign = h1_row_sign[static_cast<std::size_t>(row)];
            REAL h1_value = signed_scaled_power<REAL>(large_diff_exponents<REAL>::h1, h1_sign, h1_scale);
            left[static_cast<std::size_t>(row * k + base + 1)] = h1_value;

            REAL h2_scale = random_dyadic_unit<REAL>(&engine);
            REAL h2_sign = h2_row_sign[static_cast<std::size_t>(row)];
            REAL h2_value = signed_scaled_power<REAL>(large_diff_exponents<REAL>::h2, h2_sign, h2_scale);
            left[static_cast<std::size_t>(row * k + base + 2)] = h2_value;
        }
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            REAL p_scale = random_dyadic_unit<REAL>(&engine);
            right[static_cast<std::size_t>(base * n + col)] = p_scale;
            right[static_cast<std::size_t>((base + 3) * n + col)] = p_scale;

            REAL h1_scale = random_dyadic_unit<REAL>(&engine);
            REAL h1_sign = h1_col_sign[static_cast<std::size_t>(col)];
            REAL h1_value = h1_sign * h1_scale;
            right[static_cast<std::size_t>((base + 1) * n + col)] = h1_value;

            REAL h2_scale = random_dyadic_unit<REAL>(&engine);
            REAL h2_sign = h2_col_sign[static_cast<std::size_t>(col)];
            REAL h2_value = h2_sign * h2_scale;
            right[static_cast<std::size_t>((base + 2) * n + col)] = h2_value;
        }
    }

    for (std::ptrdiff_t pos = 4 * blocks; pos < k; ++pos) {
        for (std::ptrdiff_t row = 0; row < m; ++row) {
            REAL h2_scale = random_dyadic_unit<REAL>(&engine);
            REAL h2_sign = h2_row_sign[static_cast<std::size_t>(row)];
            REAL h2_value = signed_scaled_power<REAL>(large_diff_exponents<REAL>::h2, h2_sign, h2_scale);
            left[static_cast<std::size_t>(row * k + pos)] = h2_value;
        }
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            REAL h2_scale = random_dyadic_unit<REAL>(&engine);
            REAL h2_sign = h2_col_sign[static_cast<std::size_t>(col)];
            REAL h2_value = h2_sign * h2_scale;
            right[static_cast<std::size_t>(pos * n + col)] = h2_value;
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
    build_random_large_diff_problem<REAL>(opt.m, opt.n, opt.k, tier_seed, left, right);

    std::vector<vmplapack::Rmidrad<REAL>> result(static_cast<std::size_t>(opt.m * opt.n));
    vmplapack::VerificationStatus status =
        vmplapack::vRgemm_point<REAL>(opt.m, opt.n, opt.k, left.data(), opt.k, right.data(), opt.n, result.data(), opt.n);

    bool all_ok = true;
    bool all_covered = true;
    REAL max_radius = A::zero();
    mpfr_prec_t metric_precision = static_cast<mpfr_prec_t>(std::max(1024L, 4L * A::precision_bits() + 256L));
    mpfrxx::mpfr_class max_midpoint_error = vmplapack::oracle::zero_at(metric_precision);
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

    std::ptrdiff_t blocks = opt.k / 4;
    std::ptrdiff_t tail = opt.k - 4 * blocks;

    std::cout << std::setprecision(output_digits<REAL>()) << std::boolalpha;
    std::cout << tier << '\n';
    std::cout << "  problem = A(" << opt.m << "x" << opt.k << ") * B(" << opt.k << "x" << opt.n << "), row-major" << '\n';
    std::cout << "  construction = random dyadic-uniform [a*2^p, b*2^h1, c*2^h2, -a*2^p] blocks" << '\n';
    std::cout << "  exponents p/h1/h2 = " << large_diff_exponents<REAL>::p << "/"
              << large_diff_exponents<REAL>::h1 << "/" << large_diff_exponents<REAL>::h2 << '\n';
    std::cout << "  seed = " << static_cast<unsigned long long>(tier_seed) << '\n';
    std::cout << "  four-term blocks = " << blocks << '\n';
    std::cout << "  tail h2 terms = " << tail << '\n';
    std::cout << "  dyadic scale bits = " << dyadic_random_bits << '\n';
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
    std::cout << "  max |C.mid - oracle midpoint| = " << max_midpoint_error << '\n';
    std::cout << std::setprecision(output_digits<REAL>());
    std::cout << "  max radius = " << max_radius << '\n';
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
