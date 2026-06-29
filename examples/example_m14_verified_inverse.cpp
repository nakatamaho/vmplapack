// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

namespace {

long parse_long_env(const char* text, const char* name) {
    char* end = nullptr;
    long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0') {
        std::cerr << "invalid " << name << "=" << text << '\n';
        std::exit(EXIT_FAILURE);
    }
    return value;
}

long mpfr_precision_from_environment() {
    const char* text = std::getenv("MPFRXX_DEFAULT_PRECISION_BITS");
    if (text == nullptr) {
        return 512L;
    }

    long precision = parse_long_env(text, "MPFRXX_DEFAULT_PRECISION_BITS");
    if (precision < static_cast<long>(MPFR_PREC_MIN)) {
        std::cerr << "invalid MPFRXX_DEFAULT_PRECISION_BITS=" << text << '\n';
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
REAL power_of_two(int exponent) {
    if constexpr (std::is_same<REAL, float>::value || std::is_same<REAL, double>::value) {
        REAL value = static_cast<REAL>(std::ldexp(1.0, exponent));
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
REAL max_radius(std::ptrdiff_t n, const std::vector<vmplapack::Rmidrad<REAL>>& boxes, std::ptrdiff_t ld) {
    using A = vmplapack::Rarith<REAL>;
    REAL value = A::zero();
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            const vmplapack::Rmidrad<REAL>& box = boxes[static_cast<std::size_t>(row * ld + col)];
            if (box.rad > value) {
                value = box.rad;
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
void build_regular_matrix(std::vector<REAL>& matrix) {
    matrix.assign(static_cast<std::size_t>(9), vmplapack::Rarith<REAL>::zero());
    matrix[0] = from_int<REAL>(4);
    matrix[1] = from_int<REAL>(2);
    matrix[2] = from_int<REAL>(-1);
    matrix[3] = from_int<REAL>(1);
    matrix[4] = from_int<REAL>(3);
    matrix[5] = from_int<REAL>(2);
    matrix[6] = from_int<REAL>(0);
    matrix[7] = from_int<REAL>(-2);
    matrix[8] = from_int<REAL>(5);
}

template <class REAL>
int near_singular_gap_exponent() {
    long precision = vmplapack::Rarith<REAL>::precision_bits();
    long exponent = -precision / 4;
    if (exponent > -6L) {
        exponent = -6L;
    }
    if (exponent < -160L) {
        exponent = -160L;
    }
    return static_cast<int>(exponent);
}

template <class REAL>
void build_near_singular_matrix(std::vector<REAL>& matrix, REAL& gap) {
    matrix.assign(static_cast<std::size_t>(4), vmplapack::Rarith<REAL>::zero());
    gap = power_of_two<REAL>(near_singular_gap_exponent<REAL>());
    matrix[0] = vmplapack::Rarith<REAL>::one();
    matrix[1] = vmplapack::Rarith<REAL>::one();
    matrix[2] = vmplapack::Rarith<REAL>::one();
    REAL perturbed = vmplapack::Rarith<REAL>::one() + gap;
    matrix[3] = perturbed;
}

template <class REAL>
void show_case(const char* case_name, std::ptrdiff_t n, const std::vector<REAL>& A_data) {
    using A = vmplapack::Rarith<REAL>;

    std::vector<vmplapack::Rmidrad<REAL>> Ainv(static_cast<std::size_t>(n * n));
    vmplapack::VerificationStatus status = vmplapack::vRgeinv<REAL>(n, A_data.data(), n, Ainv.data(), n);

    std::vector<REAL> inverse_mid;
    std::vector<REAL> inverse_rad;
    collect_mid_radius(n, Ainv, n, inverse_mid, inverse_rad);

    std::cout << "  case = " << case_name << '\n';
    std::cout << "  n = " << n << '\n';
    std::cout << "  A = ";
    print_matrix(n, n, A_data, n);
    std::cout << '\n';
    std::cout << "  return status = " << status_name(status) << '\n';
    std::cout << "  all inverse boxes ok = " << std::boolalpha << all_boxes_ok(n, Ainv, n) << '\n';
    std::cout << "  inverse.mid = ";
    print_matrix(n, n, inverse_mid, n);
    std::cout << '\n';
    std::cout << "  inverse.rad = ";
    print_matrix(n, n, inverse_rad, n);
    std::cout << '\n';
    REAL radius = max_radius(n, Ainv, n);
    REAL radius_units = radius / A::unit_roundoff();
    std::cout << "  max inverse radius = " << radius << '\n';
    std::cout << std::scientific << std::setprecision(summary_digits());
    std::cout << "  max inverse radius / u = " << radius_units << '\n';
    std::cout << std::defaultfloat << std::setprecision(output_digits<REAL>());

    if (status == vmplapack::VerificationStatus::Verified) {
        std::vector<vmplapack::Rmidrad<REAL>> product(static_cast<std::size_t>(n * n));
        vmplapack::VerificationStatus product_status =
            vmplapack::vRgemm_point<REAL>(n, n, n, A_data.data(), n, inverse_mid.data(), n, product.data(), n);
        std::vector<REAL> product_mid;
        std::vector<REAL> product_rad;
        collect_mid_radius(n, product, n, product_mid, product_rad);
        REAL defect = product_identity_defect(n, product, n);
        std::cout << "  A * inverse.mid status = " << status_name(product_status) << '\n';
        std::cout << "  A * inverse.mid midpoint = ";
        print_matrix(n, n, product_mid, n);
        std::cout << '\n';
        std::cout << "  A * inverse.mid radius = ";
        print_matrix(n, n, product_rad, n);
        std::cout << '\n';
        std::cout << "  verified max |I - A*inverse.mid| bound = " << defect << '\n';
    }

    const vmplapack::Rmidrad<REAL>& first = Ainv[0];
    std::cout << "  inverse(0,0) interval = [" << lower_display(first) << ", " << upper_display(first) << "]" << '\n';
    std::cout << "  inverse(0,0) status = " << status_name(first.status) << '\n';
}

template <class REAL>
void show_tier(const char* tier) {
    std::cout << std::setprecision(output_digits<REAL>()) << std::boolalpha;
    std::cout << tier << '\n';

    std::vector<REAL> regular;
    build_regular_matrix(regular);
    show_case("regular 3x3", static_cast<std::ptrdiff_t>(3), regular);

    std::vector<REAL> near_singular;
    REAL gap = vmplapack::Rarith<REAL>::zero();
    build_near_singular_matrix(near_singular, gap);
    std::cout << "  near-singular gap = 2^(" << near_singular_gap_exponent<REAL>() << ") = " << gap << '\n';
    show_case("near-singular 2x2", static_cast<std::ptrdiff_t>(2), near_singular);
}

} // namespace

int main() {
    show_tier<float>("float");
    show_tier<double>("double");
    long mpfr_precision = mpfr_precision_from_environment();
    mpfrxx::set_default_precision_bits(static_cast<mpfr_prec_t>(mpfr_precision));
    std::string label = std::string("mpfr@") + std::to_string(mpfr_precision);
    show_tier<mpfrxx::mpfr_class>(label.c_str());
    return EXIT_SUCCESS;
}
