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
vmplapack::Rmidrad<REAL> sentinel_box() {
    vmplapack::Rmidrad<REAL> box;
    box.mid = from_int<REAL>(99);
    box.rad = from_int<REAL>(99);
    box.status = vmplapack::Rstatus::invalid_input;
    return box;
}

template <class REAL>
void test_regular_inverse_layout(const char* tier) {
    // Catches row-major/leading-dimension mistakes and certificates that miss exact inverse entries.
    std::vector<REAL> A_data(static_cast<std::size_t>(6), from_int<REAL>(77));
    A_data[0] = from_int<REAL>(2);
    A_data[1] = from_int<REAL>(1);
    A_data[3] = from_int<REAL>(1);
    A_data[4] = from_int<REAL>(1);

    std::vector<vmplapack::Rmidrad<REAL>> Ainv(static_cast<std::size_t>(6), sentinel_box<REAL>());
    vmplapack::VerificationStatus status = vmplapack::vRgeinv<REAL>(2, A_data.data(), 3, Ainv.data(), 3);

    require<REAL>(status == vmplapack::VerificationStatus::Verified, tier, "vRgeinv did not verify regular inverse");
    require<REAL>(box_covers_value(Ainv[0], from_int<REAL>(1)), tier, "Ainv(0,0) missed exact inverse");
    require<REAL>(box_covers_value(Ainv[1], from_int<REAL>(-1)), tier, "Ainv(0,1) missed exact inverse");
    require<REAL>(box_covers_value(Ainv[3], from_int<REAL>(-1)), tier, "Ainv(1,0) missed exact inverse");
    require<REAL>(box_covers_value(Ainv[4], from_int<REAL>(2)), tier, "Ainv(1,1) missed exact inverse");
    require<REAL>(Ainv[2].status == vmplapack::Rstatus::invalid_input && Ainv[5].status == vmplapack::Rstatus::invalid_input,
                  tier,
                  "vRgeinv wrote outside the active row-major inverse columns");
}

template <class REAL>
void test_boundary_rules(const char* tier) {
    // Catches drift from the frozen M14a InvalidInput and zero-dimension boundary rules.
    REAL one = vmplapack::Rarith<REAL>::one();
    vmplapack::Rmidrad<REAL> out;
    require<REAL>(vmplapack::vRgeinv<REAL>(0, nullptr, 0, nullptr, 0) == vmplapack::VerificationStatus::Verified,
                  tier,
                  "n==0 vRgeinv should write nothing and verify");
    require<REAL>(vmplapack::vRgeinv<REAL>(-1, &one, 1, &out, 1) == vmplapack::VerificationStatus::InvalidInput,
                  tier,
                  "negative n was not InvalidInput");
    require<REAL>(vmplapack::vRgeinv<REAL>(1, nullptr, 1, &out, 1) == vmplapack::VerificationStatus::InvalidInput,
                  tier,
                  "null A was not InvalidInput");
    require<REAL>(vmplapack::vRgeinv<REAL>(1, &one, 1, nullptr, 1) == vmplapack::VerificationStatus::InvalidInput,
                  tier,
                  "null Ainv was not InvalidInput");
    require<REAL>(vmplapack::vRgeinv<REAL>(1, &one, 0, &out, 1) == vmplapack::VerificationStatus::InvalidInput,
                  tier,
                  "too-small lda was not InvalidInput");
    require<REAL>(vmplapack::vRgeinv<REAL>(1, &one, 1, &out, 0) == vmplapack::VerificationStatus::InvalidInput,
                  tier,
                  "too-small ldinv was not InvalidInput");
}

template <class REAL>
void test_non_finite_inputs(const char* tier) {
    // Catches false unbounded/Verified results when non-finite A must produce non_finite boxes.
    std::vector<REAL> A_data(static_cast<std::size_t>(4), vmplapack::Rarith<REAL>::zero());
    A_data[0] = infinity_value<REAL>();
    A_data[1] = from_int<REAL>(0);
    A_data[2] = from_int<REAL>(0);
    A_data[3] = from_int<REAL>(1);
    std::vector<vmplapack::Rmidrad<REAL>> Ainv(static_cast<std::size_t>(4));

    vmplapack::VerificationStatus status = vmplapack::vRgeinv<REAL>(2, A_data.data(), 2, Ainv.data(), 2);
    require<REAL>(status == vmplapack::VerificationStatus::Unverified, tier, "non-finite vRgeinv was not Unverified");
    require<REAL>(Ainv[0].status == vmplapack::Rstatus::non_finite && Ainv[1].status == vmplapack::Rstatus::non_finite &&
                      Ainv[2].status == vmplapack::Rstatus::non_finite && Ainv[3].status == vmplapack::Rstatus::non_finite,
                  tier,
                  "non-finite vRgeinv did not write non_finite boxes");
}

template <class REAL>
void test_singular_certificate_failure(const char* tier) {
    // Catches a false inverse certificate for a singular point matrix.
    std::vector<REAL> A_data(static_cast<std::size_t>(4), vmplapack::Rarith<REAL>::zero());
    A_data[0] = from_int<REAL>(1);
    A_data[1] = from_int<REAL>(2);
    A_data[2] = from_int<REAL>(2);
    A_data[3] = from_int<REAL>(4);
    std::vector<vmplapack::Rmidrad<REAL>> Ainv(static_cast<std::size_t>(4));

    vmplapack::VerificationStatus status = vmplapack::vRgeinv<REAL>(2, A_data.data(), 2, Ainv.data(), 2);
    require<REAL>(status == vmplapack::VerificationStatus::Unverified, tier, "singular vRgeinv was not Unverified");
    require<REAL>(Ainv[0].status == vmplapack::Rstatus::unbounded && Ainv[1].status == vmplapack::Rstatus::unbounded &&
                      Ainv[2].status == vmplapack::Rstatus::unbounded && Ainv[3].status == vmplapack::Rstatus::unbounded,
                  tier,
                  "singular vRgeinv did not write unbounded boxes");
}

template <class REAL>
void test_tier(const char* tier) {
    test_regular_inverse_layout<REAL>(tier);
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
