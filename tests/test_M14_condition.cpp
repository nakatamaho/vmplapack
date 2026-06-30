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
REAL from_ratio(int numerator, int denominator) {
    REAL top = from_int<REAL>(numerator);
    REAL bottom = from_int<REAL>(denominator);
    REAL result = top / bottom;
    return result;
}

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
void test_diagonal_condition_bounds(const char* tier) {
    // Catches wrong bound direction and row-major leading-dimension padding reads.
    std::vector<REAL> A_data(static_cast<std::size_t>(6), from_int<REAL>(99));
    A_data[0] = from_int<REAL>(4);
    A_data[1] = from_int<REAL>(0);
    A_data[3] = from_int<REAL>(0);
    A_data[4] = from_int<REAL>(1);

    vmplapack::Rcond_bounds<REAL> bounds = vmplapack::vRgecon<REAL>(2, A_data.data(), 3);
    require<REAL>(bounds.status == vmplapack::Rstatus::ok, tier, "diagonal vRgecon did not return ok");
    require<REAL>(bounds.normA_upper >= from_int<REAL>(4) && bounds.normA_upper < from_int<REAL>(5),
                  tier,
                  "normA_upper is not a tight enough upper bound for diag(4,1)");
    require<REAL>(bounds.normAinv_upper >= from_int<REAL>(1) && bounds.normAinv_upper < from_int<REAL>(2),
                  tier,
                  "normAinv_upper is not a tight enough inverse-norm upper bound for diag(4,1)");
    require<REAL>(bounds.kappa_upper >= from_int<REAL>(4) && bounds.kappa_upper < from_int<REAL>(5),
                  tier,
                  "kappa_upper is not an upper bound for cond_inf(diag(4,1))");
    require<REAL>(bounds.rcond_lower <= from_ratio<REAL>(1, 4) && bounds.rcond_lower > from_ratio<REAL>(1, 5),
                  tier,
                  "rcond_lower is not a lower bound for rcond_inf(diag(4,1))");
}

template <class REAL>
void test_empty_condition_bounds(const char* tier) {
    // Catches drift from the documented empty-matrix convention for scalar condition bounds.
    vmplapack::Rcond_bounds<REAL> bounds = vmplapack::vRgecon<REAL>(0, nullptr, 0);
    require<REAL>(bounds.status == vmplapack::Rstatus::ok, tier, "n==0 vRgecon did not return ok");
    require<REAL>(bounds.normA_upper == vmplapack::Rarith<REAL>::zero() &&
                      bounds.normAinv_upper == vmplapack::Rarith<REAL>::zero(),
                  tier,
                  "n==0 vRgecon norm convention changed");
    require<REAL>(bounds.kappa_upper == vmplapack::Rarith<REAL>::one() &&
                      bounds.rcond_lower == vmplapack::Rarith<REAL>::one(),
                  tier,
                  "n==0 vRgecon condition convention changed");
}

template <class REAL>
void test_invalid_inputs(const char* tier) {
    // Catches boundary-rule drift from explicit invalid_input to undefined behavior.
    REAL one = vmplapack::Rarith<REAL>::one();
    require<REAL>(vmplapack::vRgecon<REAL>(-1, &one, 1).status == vmplapack::Rstatus::invalid_input,
                  tier,
                  "negative n was not invalid_input");
    require<REAL>(vmplapack::vRgecon<REAL>(1, nullptr, 1).status == vmplapack::Rstatus::invalid_input,
                  tier,
                  "null A was not invalid_input");
    require<REAL>(vmplapack::vRgecon<REAL>(1, &one, 0).status == vmplapack::Rstatus::invalid_input,
                  tier,
                  "too-small lda was not invalid_input");
}

template <class REAL>
void test_non_finite_input(const char* tier) {
    // Catches finite-condition claims for NaN/Inf input matrices.
    std::vector<REAL> A_data(static_cast<std::size_t>(4), vmplapack::Rarith<REAL>::zero());
    A_data[0] = infinity_value<REAL>();
    A_data[1] = from_int<REAL>(0);
    A_data[2] = from_int<REAL>(0);
    A_data[3] = from_int<REAL>(1);

    vmplapack::Rcond_bounds<REAL> bounds = vmplapack::vRgecon<REAL>(2, A_data.data(), 2);
    require<REAL>(bounds.status == vmplapack::Rstatus::non_finite, tier, "non-finite input was not non_finite");
    require<REAL>(!vmplapack::Rarith<REAL>::is_finite(bounds.kappa_upper) &&
                      bounds.rcond_lower == vmplapack::Rarith<REAL>::zero(),
                  tier,
                  "non-finite vRgecon did not return infinite/zero diagnostic bounds");
}

template <class REAL>
void test_singular_certificate_failure(const char* tier) {
    // Catches a false finite condition certificate for a singular point matrix.
    std::vector<REAL> A_data(static_cast<std::size_t>(4), vmplapack::Rarith<REAL>::zero());
    A_data[0] = from_int<REAL>(1);
    A_data[1] = from_int<REAL>(2);
    A_data[2] = from_int<REAL>(2);
    A_data[3] = from_int<REAL>(4);

    vmplapack::Rcond_bounds<REAL> bounds = vmplapack::vRgecon<REAL>(2, A_data.data(), 2);
    require<REAL>(bounds.status == vmplapack::Rstatus::unbounded, tier, "singular vRgecon was not unbounded");
    require<REAL>(bounds.normA_upper >= from_int<REAL>(6), tier, "singular vRgecon lost the finite A-norm diagnostic");
    require<REAL>(!vmplapack::Rarith<REAL>::is_finite(bounds.kappa_upper) &&
                      bounds.rcond_lower == vmplapack::Rarith<REAL>::zero(),
                  tier,
                  "singular vRgecon did not return infinite/zero diagnostic bounds");
}

template <class REAL>
void test_tier(const char* tier) {
    test_diagonal_condition_bounds<REAL>(tier);
    test_empty_condition_bounds<REAL>(tier);
    test_invalid_inputs<REAL>(tier);
    test_non_finite_input<REAL>(tier);
    test_singular_certificate_failure<REAL>(tier);
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
