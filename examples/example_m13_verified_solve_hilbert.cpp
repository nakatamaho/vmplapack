// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

namespace {

struct options {
    std::ptrdiff_t n;
    bool no_matrices;
};

options default_options() {
    options opt;
    opt.n = 10;
    opt.no_matrices = false;
    return opt;
}

void print_usage(const char* program) {
    std::cout << "usage: " << program << " [options]\n"
              << "  --n N             Hilbert size, N >= 1 (default 10)\n"
              << "  --no-matrices     suppress H, x_seed, b, x.mid, and x.rad\n"
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
            long parsed = 0;
            if (!parse_long(argv[++i], &parsed)) {
                std::cerr << "invalid --n\n";
                return false;
            }
            opt->n = static_cast<std::ptrdiff_t>(parsed);
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
    return 80;
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
REAL reciprocal_int(int denominator) {
    REAL one = vmplapack::Rarith<REAL>::one();
    REAL den = from_int<REAL>(denominator);
    REAL value = one / den;
    return value;
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
REAL abs_diff_up(REAL a, REAL b) {
    typename vmplapack::Rarith<REAL>::round_up scope;
    REAL diff = (a >= b) ? (a - b) : (b - a);
    return diff;
}

template <class REAL>
void build_hilbert_problem(std::ptrdiff_t n,
                           std::vector<REAL>& H,
                           std::vector<REAL>& x_seed,
                           std::vector<REAL>& b) {
    using A = vmplapack::Rarith<REAL>;

    H.assign(static_cast<std::size_t>(n * n), A::zero());
    x_seed.assign(static_cast<std::size_t>(n), A::zero());
    b.assign(static_cast<std::size_t>(n), A::zero());

    for (std::ptrdiff_t row = 0; row < n; ++row) {
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            int denominator = static_cast<int>(row + col + 1);
            H[static_cast<std::size_t>(row * n + col)] = reciprocal_int<REAL>(denominator);
        }
    }

    for (std::ptrdiff_t i = 0; i < n; ++i) {
        int magnitude = static_cast<int>((i % 4) + 1);
        int signed_value = (i % 2 == 0) ? magnitude : -magnitude;
        x_seed[static_cast<std::size_t>(i)] = from_int<REAL>(signed_value);
    }

    for (std::ptrdiff_t row = 0; row < n; ++row) {
        b[static_cast<std::size_t>(row)] = vmplapack::Rdot<REAL>(n, H.data() + row * n, 1, x_seed.data(), 1);
    }
}

template <class REAL>
REAL max_vector_value(const std::vector<REAL>& values) {
    using A = vmplapack::Rarith<REAL>;
    REAL max_value = A::zero();
    for (std::size_t i = 0; i < values.size(); ++i) {
        REAL abs_value = A::abs(values[i]);
        if (abs_value > max_value) {
            max_value = abs_value;
        }
    }
    return max_value;
}

template <class REAL>
void show_tier(const char* tier, const options& opt) {
    using A = vmplapack::Rarith<REAL>;

    std::vector<REAL> H;
    std::vector<REAL> x_seed;
    std::vector<REAL> b;
    build_hilbert_problem(opt.n, H, x_seed, b);

    std::vector<vmplapack::Rmidrad<REAL>> result(static_cast<std::size_t>(opt.n));
    vmplapack::VerificationStatus status =
        vmplapack::vRgesv<REAL>(opt.n, H.data(), opt.n, b.data(), 1, result.data(), 1);

    bool all_ok = true;
    REAL max_radius = A::zero();
    REAL max_mid_error = A::zero();
    std::ptrdiff_t worst_radius_index = 0;
    std::ptrdiff_t worst_error_index = 0;
    std::vector<REAL> midpoints(static_cast<std::size_t>(opt.n), A::zero());
    std::vector<REAL> radii(static_cast<std::size_t>(opt.n), A::zero());
    for (std::ptrdiff_t i = 0; i < opt.n; ++i) {
        const vmplapack::Rmidrad<REAL>& box = result[static_cast<std::size_t>(i)];
        midpoints[static_cast<std::size_t>(i)] = box.mid;
        radii[static_cast<std::size_t>(i)] = box.rad;
        if (box.status != vmplapack::Rstatus::ok) {
            all_ok = false;
        }
        if (box.rad > max_radius) {
            max_radius = box.rad;
            worst_radius_index = i;
        }
        REAL error = abs_diff_up(box.mid, x_seed[static_cast<std::size_t>(i)]);
        if (error > max_mid_error) {
            max_mid_error = error;
            worst_error_index = i;
        }
    }

    std::vector<vmplapack::Rmidrad<REAL>> residual(static_cast<std::size_t>(opt.n));
    vmplapack::Rstatus residual_status =
        vmplapack::vRresidual<REAL>(opt.n, opt.n, H.data(), opt.n, midpoints.data(), b.data(), residual.data());
    std::vector<REAL> residual_mid(static_cast<std::size_t>(opt.n), A::zero());
    std::vector<REAL> residual_rad(static_cast<std::size_t>(opt.n), A::zero());
    for (std::ptrdiff_t i = 0; i < opt.n; ++i) {
        residual_mid[static_cast<std::size_t>(i)] = residual[static_cast<std::size_t>(i)].mid;
        residual_rad[static_cast<std::size_t>(i)] = residual[static_cast<std::size_t>(i)].rad;
    }

    REAL max_abs_seed = max_vector_value(x_seed);
    REAL relative_radius = A::zero();
    if (max_abs_seed != A::zero()) {
        relative_radius = max_radius / max_abs_seed;
    }
    REAL radius_units = max_radius / A::unit_roundoff();
    REAL error_units = max_mid_error / A::unit_roundoff();

    std::cout << std::setprecision(output_digits<REAL>()) << std::boolalpha;
    std::cout << tier << '\n';
    std::cout << "  construction = Hilbert point system, b = Rdot(H, x_seed) in this tier" << '\n';
    std::cout << "  n = " << opt.n << '\n';
    if (!opt.no_matrices && opt.n <= 10) {
        std::cout << "  H = ";
        print_matrix(opt.n, opt.n, H, opt.n);
        std::cout << '\n';
        std::cout << "  x_seed = ";
        print_vector(x_seed);
        std::cout << '\n';
        std::cout << "  b = ";
        print_vector(b);
        std::cout << '\n';
        std::cout << "  x.mid = ";
        print_vector(midpoints);
        std::cout << '\n';
        std::cout << "  x.rad = ";
        print_vector(radii);
        std::cout << '\n';
    }
    std::cout << "  return status = " << status_name(status) << '\n';
    std::cout << "  all components ok = " << all_ok << '\n';
    std::cout << "  residual status at x.mid = " << status_name(residual_status) << '\n';
    std::cout << "  residual.mid = ";
    print_vector(residual_mid);
    std::cout << '\n';
    std::cout << "  residual.rad = ";
    print_vector(residual_rad);
    std::cout << '\n';
    std::cout << "  max |mid - x_seed| = " << max_mid_error << '\n';
    std::cout << "  max radius = " << max_radius << '\n';
    std::cout << "  relative radius vs max |x_seed| = " << relative_radius << '\n';
    std::cout << std::scientific << std::setprecision(summary_digits());
    std::cout << "  max |mid - x_seed| / u = " << error_units << '\n';
    std::cout << "  max radius / u = " << radius_units << '\n';
    std::cout << std::defaultfloat << std::setprecision(output_digits<REAL>());
    const vmplapack::Rmidrad<REAL>& worst_radius = result[static_cast<std::size_t>(worst_radius_index)];
    const vmplapack::Rmidrad<REAL>& worst_error = result[static_cast<std::size_t>(worst_error_index)];
    std::cout << "  worst radius component = x(" << worst_radius_index << ")" << '\n';
    std::cout << "    seed = " << x_seed[static_cast<std::size_t>(worst_radius_index)] << '\n';
    std::cout << "    mid = " << worst_radius.mid << '\n';
    std::cout << "    rad = " << worst_radius.rad << '\n';
    std::cout << "    interval = [" << lower_display(worst_radius) << ", " << upper_display(worst_radius) << "]" << '\n';
    std::cout << "    status = " << status_name(worst_radius.status) << '\n';
    std::cout << "  worst midpoint drift component = x(" << worst_error_index << ")" << '\n';
    std::cout << "    seed = " << x_seed[static_cast<std::size_t>(worst_error_index)] << '\n';
    std::cout << "    mid = " << worst_error.mid << '\n';
    std::cout << "    |mid-seed| = " << max_mid_error << '\n';
}

} // namespace

int main(int argc, char** argv) {
    options opt = default_options();
    if (!parse_options(argc, argv, &opt) || !validate_options(opt)) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    show_tier<float>("float", opt);
    show_tier<double>("double", opt);
    long mpfr_precision = mpfr_precision_from_environment();
    mpfrxx::set_default_precision_bits(static_cast<mpfr_prec_t>(mpfr_precision));
    std::string mpfr_label = std::string("mpfr@") + std::to_string(mpfr_precision);
    show_tier<mpfrxx::mpfr_class>(mpfr_label.c_str(), opt);
    return EXIT_SUCCESS;
}
