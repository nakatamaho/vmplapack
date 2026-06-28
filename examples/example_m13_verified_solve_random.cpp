// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

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
    std::uint64_t seed;
    bool no_matrices;
};

options default_options() {
    options opt;
    opt.n = 5;
    opt.seed = 41U;
    opt.no_matrices = false;
    return opt;
}

void print_usage(const char* program) {
    std::cout << "usage: " << program << " [options]\n"
              << "  --n N             system size, N >= 1 (default 5)\n"
              << "  --seed S          random seed (default 41)\n"
              << "  --no-matrices     suppress A, b, and x_true matrices\n"
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
        if (std::strcmp(argv[i], "--no-matrices") == 0) {
            opt->no_matrices = true;
            continue;
        }
        if (i + 1 >= argc) {
            std::cerr << "missing value after " << argv[i] << '\n';
            return false;
        }
        if (std::strcmp(argv[i], "--n") == 0) {
            long value = 0;
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
    return true;
}

long mpfr_precision_from_environment() {
    const char* text = std::getenv("MPFRXX_DEFAULT_PRECISION_BITS");
    if (text == nullptr) {
        return 512L;
    }

    long precision = 0L;
    if (!parse_long(text, &precision) || precision < static_cast<long>(MPFR_PREC_MIN)) {
        std::cerr << "invalid MPFRXX_DEFAULT_PRECISION_BITS=" << text << "\n";
        std::exit(EXIT_FAILURE);
    }
    return precision;
}

template <class REAL>
int output_digits() {
    return std::numeric_limits<REAL>::max_digits10;
}

template <>
int output_digits<mpfrxx::mpfr_class>() {
    long bits = vmplapack::Rarith<mpfrxx::mpfr_class>::precision_bits();
    long digits = ((bits * 30103L) + 99999L) / 100000L + 8L;
    if (digits < 20L) {
        digits = 20L;
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
REAL abs_diff_up(REAL a, REAL b) {
    typename vmplapack::Rarith<REAL>::round_up scope;
    REAL diff = (a >= b) ? (a - b) : (b - a);
    return diff;
}

template <class REAL>
void build_problem(std::ptrdiff_t n,
                   std::uint64_t seed,
                   std::vector<REAL>& A_data,
                   std::vector<REAL>& x_true,
                   std::vector<REAL>& b) {
    using A = vmplapack::Rarith<REAL>;

    std::mt19937_64 engine(seed);
    A_data.assign(static_cast<std::size_t>(n * n), A::zero());
    x_true.assign(static_cast<std::size_t>(n), A::zero());
    b.assign(static_cast<std::size_t>(n), A::zero());

    for (std::ptrdiff_t i = 0; i < n; ++i) {
        int raw = static_cast<int>(engine() % 7ULL) - 3;
        if (raw == 0) {
            raw = 1;
        }
        x_true[static_cast<std::size_t>(i)] = from_int<REAL>(raw);
    }

    for (std::ptrdiff_t row = 0; row < n; ++row) {
        int row_abs_sum = 0;
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            if (row == col) {
                continue;
            }
            int value = static_cast<int>(engine() % 7ULL) - 3;
            A_data[static_cast<std::size_t>(row * n + col)] = from_int<REAL>(value);
            row_abs_sum += (value < 0) ? -value : value;
        }
        int diagonal = row_abs_sum + static_cast<int>(n) + 3 + static_cast<int>(engine() % 3ULL);
        A_data[static_cast<std::size_t>(row * n + row)] = from_int<REAL>(diagonal);
    }

    for (std::ptrdiff_t row = 0; row < n; ++row) {
        REAL acc = A::zero();
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            REAL prod = A_data[static_cast<std::size_t>(row * n + col)] * x_true[static_cast<std::size_t>(col)];
            REAL next = acc + prod;
            acc = next;
        }
        b[static_cast<std::size_t>(row)] = acc;
    }
}

template <class REAL>
void show_tier(const char* tier, const options& opt, std::uint64_t tier_seed) {
    using A = vmplapack::Rarith<REAL>;

    std::vector<REAL> A_data;
    std::vector<REAL> x_true;
    std::vector<REAL> b;
    build_problem(opt.n, tier_seed, A_data, x_true, b);

    std::vector<vmplapack::Rmidrad<REAL>> result(static_cast<std::size_t>(opt.n));
    vmplapack::VerificationStatus status =
        vmplapack::vRgesv<REAL>(opt.n, A_data.data(), opt.n, b.data(), 1, result.data(), 1);

    bool all_ok = true;
    bool all_covered = true;
    REAL max_radius = A::zero();
    REAL max_mid_error = A::zero();
    std::ptrdiff_t worst_index = 0;
    for (std::ptrdiff_t i = 0; i < opt.n; ++i) {
        const vmplapack::Rmidrad<REAL>& box = result[static_cast<std::size_t>(i)];
        if (box.status != vmplapack::Rstatus::ok) {
            all_ok = false;
        }
        if (!box_covers_value(box, x_true[static_cast<std::size_t>(i)])) {
            all_covered = false;
        }
        if (box.rad > max_radius) {
            max_radius = box.rad;
            worst_index = i;
        }
        REAL err = abs_diff_up(box.mid, x_true[static_cast<std::size_t>(i)]);
        if (err > max_mid_error) {
            max_mid_error = err;
        }
    }

    REAL radius_units = max_radius / A::unit_roundoff();
    REAL error_units = max_mid_error / A::unit_roundoff();

    std::cout << std::setprecision(output_digits<REAL>()) << std::boolalpha;
    std::cout << tier << '\n';
    std::cout << "  construction = random integer strictly diagonally dominant system" << '\n';
    std::cout << "  n = " << opt.n << '\n';
    std::cout << "  seed = " << static_cast<unsigned long long>(tier_seed) << '\n';
    if (!opt.no_matrices && opt.n <= 8) {
        std::cout << "  A = ";
        print_matrix(opt.n, opt.n, A_data, opt.n);
        std::cout << '\n';
        std::cout << "  x_true = ";
        print_vector(x_true);
        std::cout << '\n';
        std::cout << "  b = ";
        print_vector(b);
        std::cout << '\n';
    }
    std::cout << "  return status = " << status_name(status) << '\n';
    std::cout << "  all components ok = " << all_ok << '\n';
    std::cout << "  all exact components covered = " << all_covered << '\n';
    std::cout << "  max |mid - x_true| = " << max_mid_error << '\n';
    std::cout << "  max radius = " << max_radius << '\n';
    std::cout << std::scientific << std::setprecision(summary_digits());
    std::cout << "  max |mid - x_true| / u = " << error_units << '\n';
    std::cout << "  max radius / u = " << radius_units << '\n';
    std::cout << std::defaultfloat << std::setprecision(output_digits<REAL>());
    const vmplapack::Rmidrad<REAL>& worst = result[static_cast<std::size_t>(worst_index)];
    std::cout << "  worst radius component = x(" << worst_index << ")" << '\n';
    std::cout << "    x_true = " << x_true[static_cast<std::size_t>(worst_index)] << '\n';
    std::cout << "    mid = " << worst.mid << '\n';
    std::cout << "    rad = " << worst.rad << '\n';
    std::cout << "    interval = [" << lower_display(worst) << ", " << upper_display(worst) << "]" << '\n';
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
    long mpfr_precision = mpfr_precision_from_environment();
    mpfrxx::set_default_precision_bits(static_cast<mpfr_prec_t>(mpfr_precision));
    std::string mpfr_label = std::string("mpfr@") + std::to_string(mpfr_precision);
    show_tier<mpfrxx::mpfr_class>(mpfr_label.c_str(), opt, opt.seed + 2003U);
    return EXIT_SUCCESS;
}
