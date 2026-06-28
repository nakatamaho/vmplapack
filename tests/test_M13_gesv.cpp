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
void fill_regular_matrix(std::vector<REAL>& A_data) {
    A_data.assign(static_cast<std::size_t>(4), vmplapack::Rarith<REAL>::zero());
    A_data[0] = from_int<REAL>(2);
    A_data[1] = from_int<REAL>(1);
    A_data[2] = from_int<REAL>(-1);
    A_data[3] = from_int<REAL>(4);
}

template <class REAL>
void test_vector_solve_strides(const char* tier) {
    // Catches vector-overload stride mistakes and certificates that do not include the exact solution.
    std::vector<REAL> A_data;
    fill_regular_matrix(A_data);

    std::vector<REAL> b(static_cast<std::size_t>(3), vmplapack::Rarith<REAL>::zero());
    b[0] = from_int<REAL>(10);
    b[2] = from_int<REAL>(13);

    std::vector<vmplapack::Rmidrad<REAL>> x(static_cast<std::size_t>(3));
    vmplapack::VerificationStatus status = vmplapack::vRgesv<REAL>(2, A_data.data(), 2, b.data(), 2, x.data(), 2);
    require<REAL>(status == vmplapack::VerificationStatus::Verified, tier, "strided vector vRgesv was not Verified");
    require<REAL>(box_covers_value(x[0], from_int<REAL>(3)), tier, "x[0] did not enclose the exact solution");
    require<REAL>(box_covers_value(x[2], from_int<REAL>(4)), tier, "x[1] did not enclose the exact solution");
}

template <class REAL>
void test_matrix_rhs_layout(const char* tier) {
    // Catches row-major leading-dimension mistakes for A, B, and X in the matrix-RHS overload.
    std::vector<REAL> A_data(static_cast<std::size_t>(6), from_int<REAL>(99));
    A_data[0] = from_int<REAL>(2);
    A_data[1] = from_int<REAL>(1);
    A_data[3] = from_int<REAL>(-1);
    A_data[4] = from_int<REAL>(4);

    std::vector<REAL> B_data(static_cast<std::size_t>(6), from_int<REAL>(77));
    B_data[0] = from_int<REAL>(10);
    B_data[1] = from_int<REAL>(1);
    B_data[3] = from_int<REAL>(13);
    B_data[4] = from_int<REAL>(22);

    std::vector<vmplapack::Rmidrad<REAL>> X_data(static_cast<std::size_t>(6));
    vmplapack::VerificationStatus status = vmplapack::vRgesv<REAL>(2, 2, A_data.data(), 3, B_data.data(), 3, X_data.data(), 3);
    require<REAL>(status == vmplapack::VerificationStatus::Verified, tier, "matrix RHS vRgesv was not Verified");
    require<REAL>(box_covers_value(X_data[0], from_int<REAL>(3)), tier, "X(0,0) did not enclose the exact solution");
    require<REAL>(box_covers_value(X_data[1], from_int<REAL>(-2)), tier, "X(0,1) did not enclose the exact solution");
    require<REAL>(box_covers_value(X_data[3], from_int<REAL>(4)), tier, "X(1,0) did not enclose the exact solution");
    require<REAL>(box_covers_value(X_data[4], from_int<REAL>(5)), tier, "X(1,1) did not enclose the exact solution");
}

template <class REAL>
void test_boundary_rules(const char* tier) {
    // Catches drift from the frozen M13 InvalidInput and zero-dimension boundary rules.
    REAL one = vmplapack::Rarith<REAL>::one();
    vmplapack::Rmidrad<REAL> out;
    require<REAL>(vmplapack::vRgesv<REAL>(0, nullptr, 0, nullptr, 1, nullptr, 1) == vmplapack::VerificationStatus::Verified,
                  tier,
                  "n==0 vector vRgesv should write nothing and verify");
    require<REAL>(vmplapack::vRgesv<REAL>(2, 0, nullptr, 0, nullptr, 0, nullptr, 0) == vmplapack::VerificationStatus::Verified,
                  tier,
                  "nrhs==0 matrix vRgesv should write nothing and verify");
    require<REAL>(vmplapack::vRgesv<REAL>(-1, &one, 1, &one, 1, &out, 1) == vmplapack::VerificationStatus::InvalidInput,
                  tier,
                  "negative n was not InvalidInput");
    require<REAL>(vmplapack::vRgesv<REAL>(1, nullptr, 1, &one, 1, &out, 1) == vmplapack::VerificationStatus::InvalidInput,
                  tier,
                  "null A was not InvalidInput");
    require<REAL>(vmplapack::vRgesv<REAL>(1, &one, 0, &one, 1, &out, 1) == vmplapack::VerificationStatus::InvalidInput,
                  tier,
                  "too-small lda was not InvalidInput");
    require<REAL>(vmplapack::vRgesv<REAL>(1, &one, 1, &one, 0, &out, 1) == vmplapack::VerificationStatus::InvalidInput,
                  tier,
                  "non-positive incb was not InvalidInput");
    require<REAL>(vmplapack::vRgesv<REAL>(1, &one, 1, &one, 1, &out, 0) == vmplapack::VerificationStatus::InvalidInput,
                  tier,
                  "non-positive incx was not InvalidInput");
    require<REAL>(vmplapack::vRgesv<REAL>(1, 2, &one, 1, &one, 1, &out, 2) == vmplapack::VerificationStatus::InvalidInput,
                  tier,
                  "too-small ldb was not InvalidInput");
    require<REAL>(vmplapack::vRgesv<REAL>(1, 2, &one, 1, &one, 2, &out, 1) == vmplapack::VerificationStatus::InvalidInput,
                  tier,
                  "too-small ldx was not InvalidInput");
}

template <class REAL>
void test_non_finite_inputs(const char* tier) {
    // Catches false unbounded/Verified results when non-finite inputs must produce non_finite boxes.
    std::vector<REAL> A_data;
    fill_regular_matrix(A_data);
    std::vector<REAL> b(static_cast<std::size_t>(2), vmplapack::Rarith<REAL>::zero());
    b[0] = from_int<REAL>(10);
    b[1] = infinity_value<REAL>();
    std::vector<vmplapack::Rmidrad<REAL>> x(static_cast<std::size_t>(2));

    vmplapack::VerificationStatus status = vmplapack::vRgesv<REAL>(2, A_data.data(), 2, b.data(), 1, x.data(), 1);
    require<REAL>(status == vmplapack::VerificationStatus::Unverified, tier, "non-finite vector vRgesv was not Unverified");
    require<REAL>(x[0].status == vmplapack::Rstatus::non_finite && x[1].status == vmplapack::Rstatus::non_finite,
                  tier,
                  "non-finite vector vRgesv did not write non_finite boxes");
}

template <class REAL>
void test_singular_certificate_failure(const char* tier) {
    // Catches a false nonsingularity certificate when the approximate inverse cannot certify alpha < 1.
    std::vector<REAL> A_data(static_cast<std::size_t>(4), vmplapack::Rarith<REAL>::zero());
    A_data[0] = from_int<REAL>(1);
    A_data[1] = from_int<REAL>(2);
    A_data[2] = from_int<REAL>(2);
    A_data[3] = from_int<REAL>(4);
    std::vector<REAL> b(static_cast<std::size_t>(2), vmplapack::Rarith<REAL>::zero());
    b[0] = from_int<REAL>(1);
    b[1] = from_int<REAL>(2);
    std::vector<vmplapack::Rmidrad<REAL>> x(static_cast<std::size_t>(2));

    vmplapack::VerificationStatus status = vmplapack::vRgesv<REAL>(2, A_data.data(), 2, b.data(), 1, x.data(), 1);
    require<REAL>(status == vmplapack::VerificationStatus::Unverified, tier, "singular vRgesv was not Unverified");
    require<REAL>(x[0].status == vmplapack::Rstatus::unbounded && x[1].status == vmplapack::Rstatus::unbounded,
                  tier,
                  "singular vRgesv did not write unbounded boxes");
}

template <class REAL>
void test_tier(const char* tier) {
    test_vector_solve_strides<REAL>(tier);
    test_matrix_rhs_layout<REAL>(tier);
    test_boundary_rules<REAL>(tier);
    test_non_finite_inputs<REAL>(tier);
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
