// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

namespace {

struct options {
    std::ptrdiff_t n;
    long float_cond_power;
    long double_cond_power;
    long mpfr_cond_power;
    std::uint64_t seed;
    bool print_matrices;
};

options default_options() {
    options opt;
    opt.n = 12;
    opt.float_cond_power = 12L;
    opt.double_cond_power = 32L;
    opt.mpfr_cond_power = 160L;
    opt.seed = 17U;
    opt.print_matrices = false;
    return opt;
}

void print_usage(const char* program) {
    std::cout << "usage: " << program << " [options]\n"
              << "  --n N               matrix size, N >= 1 (default 12)\n"
              << "  --seed S            random seed for Householder signs (default 17)\n"
              << "  --cond P            set all target condition powers to 2^P\n"
              << "  --float-cond P      float target condition power (default 12)\n"
              << "  --double-cond P     double target condition power (default 32)\n"
              << "  --mpfr-cond P       MPFR target condition power (default 160)\n"
              << "  --print-matrices    print A, inverse.mid, and inverse.rad\n"
              << "environment:\n"
              << "  MPFRXX_DEFAULT_PRECISION_BITS sets the MPFR tier precision (default 512)\n";
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
        if (std::strcmp(argv[i], "--print-matrices") == 0) {
            opt->print_matrices = true;
            continue;
        }
        if (i + 1 >= argc) {
            std::cerr << "missing value after " << argv[i] << '\n';
            return false;
        }
        if (std::strcmp(argv[i], "--n") == 0) {
            long value = 0L;
            if (!parse_long(argv[++i], &value)) {
                std::cerr << "invalid --n\n";
                return false;
            }
            opt->n = static_cast<std::ptrdiff_t>(value);
        } else if (std::strcmp(argv[i], "--seed") == 0) {
            std::uint64_t value = 0U;
            if (!parse_u64(argv[++i], &value)) {
                std::cerr << "invalid --seed\n";
                return false;
            }
            opt->seed = value;
        } else if (std::strcmp(argv[i], "--cond") == 0) {
            long value = 0L;
            if (!parse_long(argv[++i], &value)) {
                std::cerr << "invalid --cond\n";
                return false;
            }
            opt->float_cond_power = value;
            opt->double_cond_power = value;
            opt->mpfr_cond_power = value;
        } else if (std::strcmp(argv[i], "--float-cond") == 0) {
            if (!parse_long(argv[++i], &opt->float_cond_power)) {
                std::cerr << "invalid --float-cond\n";
                return false;
            }
        } else if (std::strcmp(argv[i], "--double-cond") == 0) {
            if (!parse_long(argv[++i], &opt->double_cond_power)) {
                std::cerr << "invalid --double-cond\n";
                return false;
            }
        } else if (std::strcmp(argv[i], "--mpfr-cond") == 0) {
            if (!parse_long(argv[++i], &opt->mpfr_cond_power)) {
                std::cerr << "invalid --mpfr-cond\n";
                return false;
            }
        } else {
            std::cerr << "unknown option: " << argv[i] << '\n';
            return false;
        }
    }
    return true;
}

bool validate_options(const options& opt) {
    if (opt.n < 1) {
        std::cerr << "requires --n >= 1\n";
        return false;
    }
    if (opt.float_cond_power < 0 || opt.double_cond_power < 0 || opt.mpfr_cond_power < 0) {
        std::cerr << "condition powers must be nonnegative\n";
        return false;
    }
    return true;
}

long mpfr_precision_from_environment() {
    const char* text = std::getenv("MPFRXX_DEFAULT_PRECISION_BITS");
    if (text == nullptr) {
        return 512L;
    }
    long value = 0L;
    if (!parse_long(text, &value) || value < static_cast<long>(MPFR_PREC_MIN)) {
        std::cerr << "invalid MPFRXX_DEFAULT_PRECISION_BITS=" << text << '\n';
        std::exit(EXIT_FAILURE);
    }
    return value;
}

template <class REAL>
int output_digits() {
    return std::numeric_limits<REAL>::max_digits10;
}

template <>
int output_digits<mpfrxx::mpfr_class>() {
    long bits = vmplapack::Rarith<mpfrxx::mpfr_class>::precision_bits();
    long digits = ((bits * 30103L) + 99999L) / 100000L + 8L;
    if (digits < 30L) {
        digits = 30L;
    }
    return static_cast<int>(digits);
}

int summary_digits() {
    return 6;
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
REAL power_of_two(long exponent) {
    if constexpr (std::is_same<REAL, float>::value || std::is_same<REAL, double>::value) {
        REAL value = static_cast<REAL>(std::ldexp(1.0, static_cast<int>(exponent)));
        return value;
    } else if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        REAL value = REAL::with_precision(static_cast<mpfr_prec_t>(vmplapack::Rarith<REAL>::precision_bits()));
        mpfr_set_ui_2exp(value.mpfr_data(), 1, static_cast<mpfr_exp_t>(exponent), MPFR_RNDN);
        return value;
    } else {
        static_assert(vmplapack::Ralways_false<REAL>::value, "power_of_two does not support this REAL type.");
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
REAL abs_diff_up(REAL a, REAL b) {
    typename vmplapack::Rarith<REAL>::round_up scope;
    REAL diff = (a >= b) ? (a - b) : (b - a);
    return diff;
}

template <class REAL>
REAL identity_component(std::ptrdiff_t row, std::ptrdiff_t col) {
    REAL value = (row == col) ? vmplapack::Rarith<REAL>::one() : vmplapack::Rarith<REAL>::zero();
    return value;
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
void make_identity(std::ptrdiff_t n, std::vector<REAL>& matrix) {
    using A = vmplapack::Rarith<REAL>;
    matrix.assign(static_cast<std::size_t>(n * n), A::zero());
    for (std::ptrdiff_t i = 0; i < n; ++i) {
        matrix[static_cast<std::size_t>(i * n + i)] = A::one();
    }
}

template <class REAL>
void apply_signed_householder_left(std::ptrdiff_t n, const std::vector<int>& signs, std::vector<REAL>& matrix) {
    using A = vmplapack::Rarith<REAL>;
    REAL two = from_int<REAL>(2);
    REAL size = from_int<REAL>(static_cast<int>(n));
    REAL tau = two / size;

    for (std::ptrdiff_t col = 0; col < n; ++col) {
        REAL dot = A::zero();
        for (std::ptrdiff_t row = 0; row < n; ++row) {
            REAL signed_value = from_int<REAL>(signs[static_cast<std::size_t>(row)]);
            REAL product = signed_value * matrix[static_cast<std::size_t>(row * n + col)];
            REAL next = dot + product;
            dot = next;
        }
        REAL scale = tau * dot;
        for (std::ptrdiff_t row = 0; row < n; ++row) {
            REAL signed_value = from_int<REAL>(signs[static_cast<std::size_t>(row)]);
            REAL correction = signed_value * scale;
            REAL next = matrix[static_cast<std::size_t>(row * n + col)] - correction;
            matrix[static_cast<std::size_t>(row * n + col)] = next;
        }
    }
}

template <class REAL>
void build_mixing_matrix(std::ptrdiff_t n, std::uint64_t seed, std::vector<REAL>& Q) {
    make_identity(n, Q);
    if (n <= 1) {
        return;
    }

    std::mt19937_64 engine(seed);
    for (int reflector = 0; reflector < 3; ++reflector) {
        std::vector<int> signs(static_cast<std::size_t>(n), 1);
        for (std::ptrdiff_t row = 0; row < n; ++row) {
            signs[static_cast<std::size_t>(row)] = ((engine() & 1ULL) == 0ULL) ? 1 : -1;
        }
        apply_signed_householder_left(n, signs, Q);
    }
}

template <class REAL>
void build_diagonal_spectrum(std::ptrdiff_t n, long cond_power, std::vector<REAL>& diagonal) {
    diagonal.assign(static_cast<std::size_t>(n), vmplapack::Rarith<REAL>::one());
    if (n <= 1) {
        return;
    }
    for (std::ptrdiff_t i = 0; i < n; ++i) {
        long exponent = -(cond_power * static_cast<long>(i)) / static_cast<long>(n - 1);
        diagonal[static_cast<std::size_t>(i)] = power_of_two<REAL>(exponent);
    }
}

template <class REAL>
void build_high_condition_matrix(std::ptrdiff_t n,
                                 long cond_power,
                                 std::uint64_t seed,
                                 std::vector<REAL>& A_data,
                                 std::vector<REAL>& diagonal) {
    using A = vmplapack::Rarith<REAL>;
    std::vector<REAL> Q;
    std::vector<REAL> QD(static_cast<std::size_t>(n * n), A::zero());
    build_mixing_matrix(n, seed, Q);
    build_diagonal_spectrum(n, cond_power, diagonal);

    for (std::ptrdiff_t row = 0; row < n; ++row) {
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            REAL scaled = Q[static_cast<std::size_t>(row * n + col)] * diagonal[static_cast<std::size_t>(col)];
            QD[static_cast<std::size_t>(row * n + col)] = scaled;
        }
    }

    A_data.assign(static_cast<std::size_t>(n * n), A::zero());
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            REAL value = vmplapack::Rdot<REAL>(n, QD.data() + row * n, 1, Q.data() + col * n, 1);
            A_data[static_cast<std::size_t>(row * n + col)] = value;
        }
    }
}

template <class REAL>
void collect_mid_radius(std::ptrdiff_t n,
                        const std::vector<vmplapack::Rmidrad<REAL>>& boxes,
                        std::ptrdiff_t ld,
                        std::vector<REAL>& mid,
                        std::vector<REAL>& rad) {
    using A = vmplapack::Rarith<REAL>;
    mid.assign(static_cast<std::size_t>(n * n), A::zero());
    rad.assign(static_cast<std::size_t>(n * n), A::zero());
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            const vmplapack::Rmidrad<REAL>& box = boxes[static_cast<std::size_t>(row * ld + col)];
            mid[static_cast<std::size_t>(row * n + col)] = box.mid;
            rad[static_cast<std::size_t>(row * n + col)] = box.rad;
        }
    }
}

template <class REAL>
REAL max_radius(std::ptrdiff_t n,
                const std::vector<vmplapack::Rmidrad<REAL>>& boxes,
                std::ptrdiff_t ld,
                std::ptrdiff_t& worst_row,
                std::ptrdiff_t& worst_col) {
    using A = vmplapack::Rarith<REAL>;
    REAL value = A::zero();
    worst_row = 0;
    worst_col = 0;
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            const vmplapack::Rmidrad<REAL>& box = boxes[static_cast<std::size_t>(row * ld + col)];
            if (box.rad > value) {
                value = box.rad;
                worst_row = row;
                worst_col = col;
            }
        }
    }
    return value;
}

template <class REAL>
bool all_boxes_ok(std::ptrdiff_t n, const std::vector<vmplapack::Rmidrad<REAL>>& boxes, std::ptrdiff_t ld) {
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            const vmplapack::Rmidrad<REAL>& box = boxes[static_cast<std::size_t>(row * ld + col)];
            if (box.status != vmplapack::Rstatus::ok) {
                return false;
            }
        }
    }
    return true;
}

template <class REAL>
REAL product_identity_defect(std::ptrdiff_t n,
                             const std::vector<vmplapack::Rmidrad<REAL>>& product,
                             std::ptrdiff_t ld) {
    using A = vmplapack::Rarith<REAL>;
    REAL defect = A::zero();
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            const vmplapack::Rmidrad<REAL>& box = product[static_cast<std::size_t>(row * ld + col)];
            REAL target = identity_component<REAL>(row, col);
            REAL diff = abs_diff_up(box.mid, target);
            REAL bound = A::zero();
            {
                typename A::round_up scope;
                bound = diff + box.rad;
            }
            if (bound > defect) {
                defect = bound;
            }
        }
    }
    return defect;
}

template <class REAL>
void show_tier(const char* tier, const options& opt, long cond_power, std::uint64_t tier_seed) {
    using A = vmplapack::Rarith<REAL>;

    std::vector<REAL> A_data;
    std::vector<REAL> diagonal;
    build_high_condition_matrix(opt.n, cond_power, tier_seed, A_data, diagonal);

    std::vector<vmplapack::Rmidrad<REAL>> inverse(static_cast<std::size_t>(opt.n * opt.n));
    vmplapack::VerificationStatus status = vmplapack::vRgeinv<REAL>(opt.n, A_data.data(), opt.n, inverse.data(), opt.n);

    std::vector<REAL> inverse_mid;
    std::vector<REAL> inverse_rad;
    collect_mid_radius(opt.n, inverse, opt.n, inverse_mid, inverse_rad);

    std::ptrdiff_t worst_row = 0;
    std::ptrdiff_t worst_col = 0;
    REAL radius = max_radius(opt.n, inverse, opt.n, worst_row, worst_col);
    REAL radius_units = radius / A::unit_roundoff();
    REAL min_diagonal = diagonal[static_cast<std::size_t>(opt.n - 1)];

    std::cout << std::setprecision(output_digits<REAL>()) << std::boolalpha;
    std::cout << tier << '\n';
    std::cout << "  construction = Q * diag(1 ... 2^-p) * Q^T, Q = 3 signed Householder reflectors\n";
    std::cout << "  n = " << opt.n << '\n';
    std::cout << "  seed = " << tier_seed << '\n';
    std::cout << "  target cond_2 = 2^" << cond_power << '\n';
    std::cout << "  smallest diagonal scale = " << min_diagonal << '\n';
    if (opt.print_matrices) {
        std::cout << "  A = ";
        print_matrix(opt.n, opt.n, A_data, opt.n);
        std::cout << '\n';
        std::cout << "  inverse.mid = ";
        print_matrix(opt.n, opt.n, inverse_mid, opt.n);
        std::cout << '\n';
        std::cout << "  inverse.rad = ";
        print_matrix(opt.n, opt.n, inverse_rad, opt.n);
        std::cout << '\n';
    }
    std::cout << "  return status = " << status_name(status) << '\n';
    std::cout << "  all inverse boxes ok = " << all_boxes_ok(opt.n, inverse, opt.n) << '\n';
    std::cout << "  max inverse radius = " << radius << '\n';
    std::cout << std::scientific << std::setprecision(summary_digits());
    std::cout << "  max inverse radius / u = " << radius_units << '\n';
    std::cout << std::defaultfloat << std::setprecision(output_digits<REAL>());

    if (status == vmplapack::VerificationStatus::Verified) {
        std::vector<vmplapack::Rmidrad<REAL>> product(static_cast<std::size_t>(opt.n * opt.n));
        vmplapack::VerificationStatus product_status =
            vmplapack::vRgemm_point<REAL>(opt.n, opt.n, opt.n, A_data.data(), opt.n, inverse_mid.data(), opt.n, product.data(), opt.n);
        REAL defect = product_identity_defect(opt.n, product, opt.n);
        std::cout << "  A * inverse.mid status = " << status_name(product_status) << '\n';
        std::cout << "  verified max |I - A*inverse.mid| bound = " << defect << '\n';
    }

    const vmplapack::Rmidrad<REAL>& worst = inverse[static_cast<std::size_t>(worst_row * opt.n + worst_col)];
    std::cout << "  worst inverse component = (" << worst_row << ", " << worst_col << ")\n";
    std::cout << "    mid = " << worst.mid << '\n';
    std::cout << "    rad = " << worst.rad << '\n';
    std::cout << "    interval = [" << lower_display(worst) << ", " << upper_display(worst) << "]\n";
    std::cout << "    status = " << status_name(worst.status) << '\n';
}

} // namespace

int main(int argc, char** argv) {
    options opt = default_options();
    if (!parse_options(argc, argv, &opt) || !validate_options(opt)) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    show_tier<float>("float", opt, opt.float_cond_power, opt.seed);
    show_tier<double>("double", opt, opt.double_cond_power, opt.seed + 1009U);
    long mpfr_precision = mpfr_precision_from_environment();
    mpfrxx::set_default_precision_bits(static_cast<mpfr_prec_t>(mpfr_precision));
    std::string label = std::string("mpfr@") + std::to_string(mpfr_precision);
    show_tier<mpfrxx::mpfr_class>(label.c_str(), opt, opt.mpfr_cond_power, opt.seed + 2003U);
    return EXIT_SUCCESS;
}
