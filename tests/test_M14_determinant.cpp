// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

namespace {

template <class REAL>
void require(bool condition, const char* tier, const char* message) {
    if (!condition) {
        std::cerr << tier << ": " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

#ifdef VMPLAPACK_ENABLE_MPFR
long parse_expected_mpfr_precision() {
    const char* text = std::getenv("MPFRXX_DEFAULT_PRECISION_BITS");
    if (text == nullptr || text[0] == '\0') {
        return 0;
    }
    char* end = nullptr;
    long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0) {
        std::cerr << "invalid MPFRXX_DEFAULT_PRECISION_BITS=" << text << '\n';
        std::exit(EXIT_FAILURE);
    }
    return value;
}
#endif

template <class REAL>
REAL from_int(int value) {
    REAL result = static_cast<REAL>(value);
    return result;
}

#ifdef VMPLAPACK_ENABLE_MPFR
template <>
mpfrxx::mpfr_class from_int<mpfrxx::mpfr_class>(int value) {
    mpfrxx::mpfr_class result =
        mpfrxx::mpfr_class::with_precision(static_cast<mpfr_prec_t>(vmplapack::Rarith<mpfrxx::mpfr_class>::precision_bits()));
    mpfr_set_si(result.mpfr_data(), static_cast<long>(value), MPFR_RNDN);
    return result;
}
#endif

template <class REAL>
REAL infinity_value() {
    if constexpr (std::is_same<REAL, float>::value || std::is_same<REAL, double>::value) {
        return std::numeric_limits<REAL>::infinity();
#ifdef VMPLAPACK_ENABLE_MPFR
    } else if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        REAL value = REAL::with_precision(static_cast<mpfr_prec_t>(vmplapack::Rarith<REAL>::precision_bits()));
        mpfr_set_inf(value.mpfr_data(), 1);
        return value;
#endif
    } else {
        static_assert(vmplapack::Ralways_false<REAL>::value, "infinity_value does not support this REAL type.");
    }
}

template <class REAL>
bool box_covers_value(const vmplapack::Rmidrad<REAL>& box, REAL value) {
    using A = vmplapack::Rarith<REAL>;
    if (box.status != vmplapack::Rstatus::ok || box.rad < A::zero()) {
        return false;
    }
    REAL lo = A::zero();
    {
        typename A::round_down scope;
        lo = box.mid - box.rad;
    }
    REAL hi = A::zero();
    {
        typename A::round_up scope;
        hi = box.mid + box.rad;
    }
    return lo <= value && value <= hi;
}

template <class REAL>
void test_2x2_leading_dimension(const char* tier) {
    // Catches row-major/lda mistakes and wrong signs in the 2x2 determinant expansion.
    std::vector<REAL> A_data(static_cast<std::size_t>(6), from_int<REAL>(99));
    A_data[0] = from_int<REAL>(2);
    A_data[1] = from_int<REAL>(1);
    A_data[3] = from_int<REAL>(3);
    A_data[4] = from_int<REAL>(4);

    vmplapack::Rmidrad<REAL> det = vmplapack::vRgedet<REAL>(2, A_data.data(), 3);
    require<REAL>(box_covers_value(det, from_int<REAL>(5)), tier, "2x2 determinant box missed exact value 5");
    require<REAL>(det.rad < from_int<REAL>(1), tier, "2x2 determinant box is unexpectedly wide");
}

template <class REAL>
void test_3x3_signed_permutations(const char* tier) {
    // Catches parity/sign errors across all six 3x3 permutations.
    std::vector<REAL> A_data(static_cast<std::size_t>(12), from_int<REAL>(77));
    A_data[0] = from_int<REAL>(1);
    A_data[1] = from_int<REAL>(2);
    A_data[2] = from_int<REAL>(3);
    A_data[4] = from_int<REAL>(0);
    A_data[5] = from_int<REAL>(4);
    A_data[6] = from_int<REAL>(5);
    A_data[8] = from_int<REAL>(1);
    A_data[9] = from_int<REAL>(0);
    A_data[10] = from_int<REAL>(6);

    vmplapack::Rmidrad<REAL> det = vmplapack::vRgedet<REAL>(3, A_data.data(), 4);
    require<REAL>(box_covers_value(det, from_int<REAL>(22)), tier, "3x3 determinant box missed exact value 22");
    require<REAL>(det.rad < from_int<REAL>(1), tier, "3x3 determinant box is unexpectedly wide");
}

template <class REAL>
void test_singular_determinant_zero(const char* tier) {
    // Catches treating a singular finite matrix as an inverse-style certificate failure.
    std::vector<REAL> A_data(static_cast<std::size_t>(4), vmplapack::Rarith<REAL>::zero());
    A_data[0] = from_int<REAL>(1);
    A_data[1] = from_int<REAL>(2);
    A_data[2] = from_int<REAL>(2);
    A_data[3] = from_int<REAL>(4);

    vmplapack::Rmidrad<REAL> det = vmplapack::vRgedet<REAL>(2, A_data.data(), 2);
    require<REAL>(box_covers_value(det, vmplapack::Rarith<REAL>::zero()), tier, "singular determinant box missed zero");
}

template <class REAL>
void test_boundary_rules(const char* tier) {
    // Catches drift from the documented scalar determinant boundary rules.
    REAL one = vmplapack::Rarith<REAL>::one();
    vmplapack::Rmidrad<REAL> empty = vmplapack::vRgedet<REAL>(0, nullptr, 0);
    require<REAL>(empty.status == vmplapack::Rstatus::ok && empty.mid == one && empty.rad == vmplapack::Rarith<REAL>::zero(),
                  tier,
                  "n==0 determinant convention changed");
    require<REAL>(vmplapack::vRgedet<REAL>(-1, &one, 1).status == vmplapack::Rstatus::invalid_input,
                  tier,
                  "negative n was not invalid_input");
    require<REAL>(vmplapack::vRgedet<REAL>(1, nullptr, 1).status == vmplapack::Rstatus::invalid_input,
                  tier,
                  "null A was not invalid_input");
    require<REAL>(vmplapack::vRgedet<REAL>(1, &one, 0).status == vmplapack::Rstatus::invalid_input,
                  tier,
                  "too-small lda was not invalid_input");
}

template <class REAL>
void test_non_finite_input(const char* tier) {
    // Catches false finite determinant certificates for NaN/Inf input matrices.
    std::vector<REAL> A_data(static_cast<std::size_t>(4), vmplapack::Rarith<REAL>::zero());
    A_data[0] = infinity_value<REAL>();
    A_data[1] = from_int<REAL>(0);
    A_data[2] = from_int<REAL>(0);
    A_data[3] = from_int<REAL>(1);

    vmplapack::Rmidrad<REAL> det = vmplapack::vRgedet<REAL>(2, A_data.data(), 2);
    require<REAL>(det.status == vmplapack::Rstatus::non_finite, tier, "non-finite input was not non_finite");
}

template <class REAL>
void test_unsupported_large_order(const char* tier) {
    // Catches accidental factorial-time work beyond the M14c reference certificate limit.
    std::ptrdiff_t n = 9;
    std::vector<REAL> A_data(static_cast<std::size_t>(n * n), vmplapack::Rarith<REAL>::zero());
    for (std::ptrdiff_t i = 0; i < n; ++i) {
        A_data[static_cast<std::size_t>(i * n + i)] = vmplapack::Rarith<REAL>::one();
    }
    vmplapack::Rmidrad<REAL> det = vmplapack::vRgedet<REAL>(n, A_data.data(), n);
    require<REAL>(det.status == vmplapack::Rstatus::unbounded, tier, "n>8 determinant did not return unbounded");
}

template <class REAL>
void test_tier(const char* tier) {
    test_2x2_leading_dimension<REAL>(tier);
    test_3x3_signed_permutations<REAL>(tier);
    test_singular_determinant_zero<REAL>(tier);
    test_boundary_rules<REAL>(tier);
    test_non_finite_input<REAL>(tier);
    test_unsupported_large_order<REAL>(tier);
}

} // namespace

int main() {
#ifdef VMPLAPACK_ENABLE_MPFR
    long expected_mpfr_precision = parse_expected_mpfr_precision();
    if (expected_mpfr_precision != 0) {
        mpfrxx::set_default_precision_bits(static_cast<mpfr_prec_t>(expected_mpfr_precision));
    }
#endif

    test_tier<float>("float");
    test_tier<double>("double");

#ifdef VMPLAPACK_ENABLE_MPFR
    if (expected_mpfr_precision == 53 || expected_mpfr_precision == 512) {
        test_tier<mpfrxx::mpfr_class>("mpfr");
    }
#endif

    return EXIT_SUCCESS;
}
