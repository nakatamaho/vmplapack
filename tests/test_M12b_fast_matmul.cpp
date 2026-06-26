// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include "Rdot_oracle.h"

#include <vmplapack/vmplapack.h>

#include <cmath>
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
    require<REAL>(midrad_covers_oracle(box, ref), tier, message);
}

template <class REAL>
void test_fast_gemv_is_tighter_than_directed(const char* tier, int exponent) {
    // Catches an M12b matvec path that silently keeps the slow directed enclosure for every component.
    using A = vmplapack::Rarith<REAL>;

    REAL big = power_of_two_value<REAL>(exponent);
    std::vector<REAL> matrix;
    matrix.push_back(big);
    matrix.push_back(A::one());
    matrix.push_back(big);

    std::vector<REAL> x;
    x.push_back(A::one());
    x.push_back(A::one());
    x.push_back(-A::one());

    std::vector<vmplapack::Rmidrad<REAL>> out(static_cast<std::size_t>(1));
    vmplapack::VerificationStatus status =
        vmplapack::vRgemv_point<REAL>(1, 3, matrix.data(), 3, x.data(), 1, out.data());
    require<REAL>(status == vmplapack::VerificationStatus::Verified, tier, "fast vRgemv_point was not verified");
    require<REAL>(out[0].status == vmplapack::Rstatus::ok, tier, "fast vRgemv_point component was not ok");

    vmplapack::Rmidrad<REAL> directed = vmplapack::vRdot<REAL>(3, matrix.data(), 1, x.data(), 1);
    require<REAL>(directed.status == vmplapack::Rstatus::ok, tier, "directed reference vRdot was not ok");
    require<REAL>(out[0].rad < directed.rad, tier, "fast vRgemv_point radius was not tighter than directed reference");
    require_box_covers_dot<REAL>(out[0], 3, matrix.data(), 1, x.data(), 1, tier,
                                 "fast vRgemv_point component did not cover oracle dot interval");
}

template <class REAL>
void test_fast_gemm_is_tighter_than_directed(const char* tier, int exponent) {
    // Catches an M12b matmul path that computes a nearest midpoint but forgets the a-priori radius.
    using A = vmplapack::Rarith<REAL>;

    REAL big = power_of_two_value<REAL>(exponent);
    std::vector<REAL> A_data;
    A_data.push_back(big);
    A_data.push_back(A::one());
    A_data.push_back(big);

    std::vector<REAL> B_data;
    B_data.push_back(A::one());
    B_data.push_back(A::one());
    B_data.push_back(-A::one());

    std::vector<vmplapack::Rmidrad<REAL>> C_data(static_cast<std::size_t>(1));
    vmplapack::VerificationStatus status =
        vmplapack::vRgemm_point<REAL>(1, 1, 3, A_data.data(), 3, B_data.data(), 1, C_data.data(), 1);
    require<REAL>(status == vmplapack::VerificationStatus::Verified, tier, "fast vRgemm_point was not verified");
    require<REAL>(C_data[0].status == vmplapack::Rstatus::ok, tier, "fast vRgemm_point component was not ok");

    vmplapack::Rmidrad<REAL> directed = vmplapack::vRdot<REAL>(3, A_data.data(), 1, B_data.data(), 1);
    require<REAL>(directed.status == vmplapack::Rstatus::ok, tier, "directed reference vRdot was not ok");
    require<REAL>(C_data[0].rad < directed.rad, tier, "fast vRgemm_point radius was not tighter than directed reference");
    require_box_covers_dot<REAL>(C_data[0], 3, A_data.data(), 1, B_data.data(), 1, tier,
                                 "fast vRgemm_point component did not cover oracle dot interval");
}

template <class REAL>
void test_native_apriori_overflow_falls_back_to_reference(const char* tier) {
    // Catches an M12b fast path that returns unbounded when the a-priori S_up bound overflows but
    // the M12a directed reference can still certify cancellation exactly.
    using A = vmplapack::Rarith<REAL>;

    REAL huge = std::numeric_limits<REAL>::max();
    std::vector<REAL> A_data;
    A_data.push_back(huge);
    A_data.push_back(huge);

    std::vector<REAL> B_data;
    B_data.push_back(A::one());
    B_data.push_back(-A::one());

    std::vector<vmplapack::Rmidrad<REAL>> C_data(static_cast<std::size_t>(1));
    vmplapack::VerificationStatus status =
        vmplapack::vRgemm_point<REAL>(1, 1, 2, A_data.data(), 2, B_data.data(), 1, C_data.data(), 1);
    require<REAL>(status == vmplapack::VerificationStatus::Verified, tier, "fallback vRgemm_point was not verified");
    require<REAL>(C_data[0].status == vmplapack::Rstatus::ok, tier, "fallback component was not ok");
    require<REAL>(C_data[0].mid == A::zero() && C_data[0].rad == A::zero(), tier,
                  "fallback component did not recover the exact zero directed reference box");
}

template <class REAL>
void run_tier(const char* tier, int exponent) {
    test_fast_gemv_is_tighter_than_directed<REAL>(tier, exponent);
    test_fast_gemm_is_tighter_than_directed<REAL>(tier, exponent);
}

} // namespace

int main() {
    run_tier<float>("float", 30);
    test_native_apriori_overflow_falls_back_to_reference<float>("float");

    run_tier<double>("double", 70);
    test_native_apriori_overflow_falls_back_to_reference<double>("double");

#ifdef VMPLAPACK_ENABLE_MPFR
    long expected_precision = parse_expected_mpfr_precision();
    if (expected_precision == 53 || expected_precision == 512) {
        mpfrxx::set_default_precision_bits(static_cast<mpfr_prec_t>(expected_precision));
        int exponent = (expected_precision == 53) ? 70 : 700;
        run_tier<mpfrxx::mpfr_class>("mpfr", exponent);
    }
#endif

    return EXIT_SUCCESS;
}
