// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include "Rdot_oracle.h"
#include "Rgendot.h"

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

template <class REAL>
REAL quiet_nan_value() {
    if constexpr (std::is_same<REAL, float>::value || std::is_same<REAL, double>::value) {
        return std::numeric_limits<REAL>::quiet_NaN();
    } else if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        REAL value = REAL::with_precision(static_cast<mpfr_prec_t>(vmplapack::Rarith<REAL>::precision_bits()));
        mpfr_set_nan(value.mpfr_data());
        return value;
    } else {
        static_assert(vmplapack::Ralways_false<REAL>::value, "quiet_nan_value does not support this REAL type.");
    }
}

template <class REAL>
REAL overflow_large_value() {
    if constexpr (std::is_same<REAL, float>::value || std::is_same<REAL, double>::value) {
        return std::numeric_limits<REAL>::max();
    } else if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        REAL value = REAL::with_precision(static_cast<mpfr_prec_t>(vmplapack::Rarith<REAL>::precision_bits()));
        mpfr_exp_t exponent = mpfr_get_emax() - 1;
        mpfr_set_ui_2exp(value.mpfr_data(), 1, exponent, MPFR_RNDN);
        return value;
    } else {
        static_assert(vmplapack::Ralways_false<REAL>::value, "overflow_large_value does not support this REAL type.");
    }
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
void require_box_covers_dot(const vmplapack::Rmidrad<REAL>& box,
                            std::ptrdiff_t n,
                            const REAL* x,
                            std::ptrdiff_t incx,
                            const REAL* y,
                            std::ptrdiff_t incy,
                            const char* tier,
                            const char* message) {
    vmplapack::oracle::Rdot_interval ref = vmplapack::oracle::Rdot_oracle(n, x, incx, y, incy);
    require<REAL>(box.status == vmplapack::Rstatus::ok, tier, message);
    require<REAL>(midrad_covers_oracle(box, ref), tier, message);
}

template <class REAL>
void fill_reference_problem(int exponent,
                            std::vector<REAL>& A_data,
                            std::vector<REAL>& x_data,
                            std::vector<REAL>& B_data) {
    using A = vmplapack::Rarith<REAL>;

    REAL one = A::one();
    REAL two = one + one;
    REAL large = vmplapack::gendot::power_of_two<REAL>(exponent);

    A_data.clear();
    A_data.push_back(large);
    A_data.push_back(one);
    A_data.push_back(-large);
    A_data.push_back(one);
    A_data.push_back(-two);
    A_data.push_back(one);

    x_data.clear();
    x_data.push_back(one);
    x_data.push_back(one);
    x_data.push_back(one);

    B_data.clear();
    B_data.push_back(one);
    B_data.push_back(one);
    B_data.push_back(one);
    B_data.push_back(-one);
    B_data.push_back(one);
    B_data.push_back(one);
}

template <class REAL>
void test_gemv_reference(const char* tier, int exponent) {
    // Catches M12a matvec implementations that rely on BLAS or nearest-rounded dot products instead
    // of directed lower/upper component enclosures.
    std::vector<REAL> A_data;
    std::vector<REAL> x_data;
    std::vector<REAL> B_data;
    fill_reference_problem(exponent, A_data, x_data, B_data);

    std::vector<vmplapack::Rmidrad<REAL>> out(static_cast<std::size_t>(2));
    vmplapack::VerificationStatus status =
        vmplapack::vRgemv_point<REAL>(2, 3, A_data.data(), 3, x_data.data(), 1, out.data());
    require<REAL>(status == vmplapack::VerificationStatus::Verified, tier, "vRgemv_point did not verify regular case");
    for (std::ptrdiff_t row = 0; row < 2; ++row) {
        require_box_covers_dot(out[static_cast<std::size_t>(row)],
                               3,
                               A_data.data() + row * 3,
                               1,
                               x_data.data(),
                               1,
                               tier,
                               "vRgemv_point component did not cover oracle dot interval");
    }
}

template <class REAL>
void test_gemv_strided_vector(const char* tier, int exponent) {
    // Catches a matvec implementation that validates incx but indexes x as if incx were always one.
    using A = vmplapack::Rarith<REAL>;

    std::vector<REAL> A_data;
    std::vector<REAL> compact_x;
    std::vector<REAL> B_data;
    fill_reference_problem(exponent, A_data, compact_x, B_data);

    std::vector<REAL> strided_x(static_cast<std::size_t>(5), A::zero());
    strided_x[0] = compact_x[0];
    strided_x[2] = compact_x[1];
    strided_x[4] = compact_x[2];

    std::vector<vmplapack::Rmidrad<REAL>> out(static_cast<std::size_t>(2));
    vmplapack::VerificationStatus status =
        vmplapack::vRgemv_point<REAL>(2, 3, A_data.data(), 3, strided_x.data(), 2, out.data());
    require<REAL>(status == vmplapack::VerificationStatus::Verified, tier, "strided vRgemv_point did not verify");
    for (std::ptrdiff_t row = 0; row < 2; ++row) {
        require_box_covers_dot(out[static_cast<std::size_t>(row)],
                               3,
                               A_data.data() + row * 3,
                               1,
                               strided_x.data(),
                               2,
                               tier,
                               "strided vRgemv_point component did not cover oracle dot interval");
    }
}

template <class REAL>
void test_gemm_reference(const char* tier, int exponent) {
    // Catches row/column stride mistakes in the M12a matrix product reference path.
    std::vector<REAL> A_data;
    std::vector<REAL> x_data;
    std::vector<REAL> B_data;
    fill_reference_problem(exponent, A_data, x_data, B_data);

    std::vector<vmplapack::Rmidrad<REAL>> C_data(static_cast<std::size_t>(4));
    vmplapack::VerificationStatus status =
        vmplapack::vRgemm_point<REAL>(2, 2, 3, A_data.data(), 3, B_data.data(), 2, C_data.data(), 2);
    require<REAL>(status == vmplapack::VerificationStatus::Verified, tier, "vRgemm_point did not verify regular case");
    for (std::ptrdiff_t row = 0; row < 2; ++row) {
        for (std::ptrdiff_t col = 0; col < 2; ++col) {
            require_box_covers_dot(C_data[static_cast<std::size_t>(row * 2 + col)],
                                   3,
                                   A_data.data() + row * 3,
                                   1,
                                   B_data.data() + col,
                                   2,
                                   tier,
                                   "vRgemm_point component did not cover oracle dot interval");
        }
    }
}

template <class REAL>
void test_zero_dimensions(const char* tier) {
    // Catches zero-length products that still dereference A/B/x or leave nonzero output boxes.
    using A = vmplapack::Rarith<REAL>;

    std::vector<vmplapack::Rmidrad<REAL>> gemv_out(static_cast<std::size_t>(2));
    vmplapack::VerificationStatus gemv_status =
        vmplapack::vRgemv_point<REAL>(2, 0, nullptr, 0, nullptr, 1, gemv_out.data());
    require<REAL>(gemv_status == vmplapack::VerificationStatus::Verified, tier, "n==0 vRgemv_point was not verified");
    require<REAL>(gemv_out[0].mid == A::zero() && gemv_out[0].rad == A::zero() && gemv_out[0].status == vmplapack::Rstatus::ok,
                  tier,
                  "n==0 vRgemv_point output[0] was not exact zero");
    require<REAL>(gemv_out[1].mid == A::zero() && gemv_out[1].rad == A::zero() && gemv_out[1].status == vmplapack::Rstatus::ok,
                  tier,
                  "n==0 vRgemv_point output[1] was not exact zero");

    std::vector<vmplapack::Rmidrad<REAL>> gemm_out(static_cast<std::size_t>(6));
    vmplapack::VerificationStatus gemm_status =
        vmplapack::vRgemm_point<REAL>(2, 3, 0, nullptr, 0, nullptr, 0, gemm_out.data(), 3);
    require<REAL>(gemm_status == vmplapack::VerificationStatus::Verified, tier, "k==0 vRgemm_point was not verified");
    for (std::size_t i = 0; i < gemm_out.size(); ++i) {
        require<REAL>(gemm_out[i].mid == A::zero() && gemm_out[i].rad == A::zero() && gemm_out[i].status == vmplapack::Rstatus::ok,
                      tier,
                      "k==0 vRgemm_point output was not exact zero");
    }

    require<REAL>(vmplapack::vRgemv_point<REAL>(0, 3, nullptr, 0, nullptr, 1, nullptr) ==
                      vmplapack::VerificationStatus::Verified,
                  tier,
                  "m==0 vRgemv_point should not require pointers");
    require<REAL>(vmplapack::vRgemm_point<REAL>(0, 3, 4, nullptr, 0, nullptr, 0, nullptr, 0) ==
                      vmplapack::VerificationStatus::Verified,
                  tier,
                  "m==0 vRgemm_point should not require pointers");
}

template <class REAL>
void test_invalid_inputs(const char* tier) {
    // Catches boundary-rule drift from explicit InvalidInput into scalar Rstatus or undefined behavior.
    using A = vmplapack::Rarith<REAL>;

    REAL one = A::one();
    vmplapack::Rmidrad<REAL> out;
    require<REAL>(vmplapack::vRgemv_point<REAL>(-1, 1, &one, 1, &one, 1, &out) ==
                      vmplapack::VerificationStatus::InvalidInput,
                  tier,
                  "negative m was not InvalidInput");
    require<REAL>(vmplapack::vRgemv_point<REAL>(1, 1, nullptr, 1, &one, 1, &out) ==
                      vmplapack::VerificationStatus::InvalidInput,
                  tier,
                  "null A was not InvalidInput");
    require<REAL>(vmplapack::vRgemv_point<REAL>(1, 1, &one, 1, &one, 0, &out) ==
                      vmplapack::VerificationStatus::InvalidInput,
                  tier,
                  "non-positive incx was not InvalidInput");
    require<REAL>(vmplapack::vRgemm_point<REAL>(1, 1, 1, &one, 0, &one, 1, &out, 1) ==
                      vmplapack::VerificationStatus::InvalidInput,
                  tier,
                  "too-small lda was not InvalidInput");
    require<REAL>(vmplapack::vRgemm_point<REAL>(1, 1, 1, &one, 1, &one, 1, nullptr, 1) ==
                      vmplapack::VerificationStatus::InvalidInput,
                  tier,
                  "null C was not InvalidInput");
}

template <class REAL>
void test_unverified_cases(const char* tier) {
    // Catches false Verified returns for non-finite inputs and finite-input overflow.
    using A = vmplapack::Rarith<REAL>;

    REAL one = A::one();
    REAL two = one + one;
    REAL nan_value = quiet_nan_value<REAL>();
    vmplapack::Rmidrad<REAL> out;
    vmplapack::VerificationStatus non_finite_status =
        vmplapack::vRgemv_point<REAL>(1, 1, &nan_value, 1, &one, 1, &out);
    require<REAL>(non_finite_status == vmplapack::VerificationStatus::Unverified,
                  tier,
                  "non-finite vRgemv_point input was not Unverified");
    require<REAL>(out.status == vmplapack::Rstatus::non_finite,
                  tier,
                  "non-finite vRgemv_point component did not report non_finite");

    REAL large = overflow_large_value<REAL>();
    vmplapack::VerificationStatus overflow_status =
        vmplapack::vRgemv_point<REAL>(1, 1, &large, 1, &two, 1, &out);
    require<REAL>(overflow_status == vmplapack::VerificationStatus::Unverified,
                  tier,
                  "overflow vRgemv_point input was not Unverified");
    require<REAL>(out.status == vmplapack::Rstatus::unbounded && !A::is_finite(out.rad),
                  tier,
                  "overflow vRgemv_point component did not report unbounded");
}

template <class REAL>
void test_tier(const char* tier, int exponent) {
    test_gemv_reference<REAL>(tier, exponent);
    test_gemv_strided_vector<REAL>(tier, exponent);
    test_gemm_reference<REAL>(tier, exponent);
    test_zero_dimensions<REAL>(tier);
    test_invalid_inputs<REAL>(tier);
    test_unverified_cases<REAL>(tier);
}

} // namespace

int main() {
    long expected_mpfr_precision = parse_expected_mpfr_precision();
    if (expected_mpfr_precision != 0) {
        mpfrxx::set_default_precision_bits(static_cast<mpfr_prec_t>(expected_mpfr_precision));
    }

    test_tier<float>("float", 30);
    test_tier<double>("double", 60);

    if (expected_mpfr_precision == 53) {
        test_tier<mpfrxx::mpfr_class>("mpfr", 60);
    } else if (expected_mpfr_precision == 512) {
        test_tier<mpfrxx::mpfr_class>("mpfr", 700);
    } else if (expected_mpfr_precision != 0) {
        std::cerr << "unsupported MPFRXX_DEFAULT_PRECISION_BITS=" << expected_mpfr_precision << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
