// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <type_traits>

#if (defined(__i386__) || defined(__x86_64__)) && defined(__SSE__)
#include <xmmintrin.h>
#endif

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif
#pragma STDC FENV_ACCESS ON
#pragma STDC FP_CONTRACT OFF
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace {

#if defined(__GNUC__) || defined(__clang__)
#define VMPLAPACK_M9_NOINLINE __attribute__((noinline))
#else
#define VMPLAPACK_M9_NOINLINE
#endif

template <class REAL>
void require(bool condition, const char* tier, const char* message) {
    if (!condition) {
        std::cerr << tier << ": " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

template <class REAL>
VMPLAPACK_M9_NOINLINE REAL runtime_add(REAL a, REAL b) {
    volatile REAL va = a;
    volatile REAL vb = b;
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : "+m"(va), "+m"(vb) :: "memory");
#endif
    REAL x = va;
    REAL y = vb;
    volatile REAL z = x + y;
    return z;
}

template <class REAL>
VMPLAPACK_M9_NOINLINE REAL runtime_mul(REAL a, REAL b) {
    volatile REAL va = a;
    volatile REAL vb = b;
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : "+m"(va), "+m"(vb) :: "memory");
#endif
    REAL x = va;
    REAL y = vb;
    volatile REAL z = x * y;
    return z;
}

template <class REAL>
REAL tier_power_of_two(long exponent) {
    using A = vmplapack::Rarith<REAL>;
    return std::ldexp(A::one(), static_cast<int>(exponent));
}

#ifdef VMPLAPACK_ENABLE_MPFR
template <>
mpfrxx::mpfr_class tier_power_of_two<mpfrxx::mpfr_class>(long exponent) {
    mpfr_prec_t precision = static_cast<mpfr_prec_t>(vmplapack::Rarith<mpfrxx::mpfr_class>::precision_bits());
    mpfrxx::mpfr_class value = mpfrxx::mpfr_class::with_precision(precision);
    mpfr_set_ui_2exp(value.mpfr_data(), 1, static_cast<mpfr_exp_t>(exponent), MPFR_RNDN);
    return value;
}
#endif

template <class REAL>
bool is_negative_zero(REAL value) {
    return value == vmplapack::Rarith<REAL>::zero() && std::signbit(value);
}

#ifdef VMPLAPACK_ENABLE_MPFR
template <>
bool is_negative_zero<mpfrxx::mpfr_class>(mpfrxx::mpfr_class value) {
    return value == vmplapack::Rarith<mpfrxx::mpfr_class>::zero() && mpfr_signbit(value.mpfr_data()) != 0;
}
#endif

template <class REAL>
void require_zero_value(REAL value, const char* tier, const char* message) {
    require<REAL>(value == vmplapack::Rarith<REAL>::zero(), tier, message);
}

template <class REAL>
REAL next_after_one_by_unit_roundoff() {
    using A = vmplapack::Rarith<REAL>;

    REAL two = A::one() + A::one();
    REAL step = two * A::unit_roundoff();
    REAL next = A::one() + step;
    return next;
}

template <class REAL>
void test_eft_normal_case(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches loss of round-to-nearest or a broken TwoSum residual on a normal half-ulp tie.
    REAL s = A::zero();
    REAL e = A::zero();
    vmplapack::Rtwosum(A::one(), A::unit_roundoff(), s, e);
    require<REAL>(s == A::one(), tier, "Rtwosum normal high part is wrong");
    require<REAL>(e == A::unit_roundoff(), tier, "Rtwosum normal residual is wrong");

    // Catches FastTwoSum violations under its abs(a) >= abs(b) precondition.
    vmplapack::Rfasttwosum(A::one(), A::unit_roundoff(), s, e);
    require<REAL>(s == A::one(), tier, "Rfasttwosum normal high part is wrong");
    require<REAL>(e == A::unit_roundoff(), tier, "Rfasttwosum normal residual is wrong");

    // Catches FMA-residual or contraction mistakes in TwoProduct on a known normal product.
    REAL next = next_after_one_by_unit_roundoff<REAL>();
    REAL eps = next - A::one();
    REAL p = A::zero();
    REAL r = A::zero();
    vmplapack::Rtwoproduct(next, next, p, r);
    REAL eps2 = eps + eps;
    REAL expected_p = A::one() + eps2;
    REAL expected_r = eps * eps;
    require<REAL>(p == expected_p, tier, "Rtwoproduct normal high part is wrong");
    require<REAL>(r == expected_r, tier, "Rtwoproduct normal residual is wrong");
}

template <class REAL>
void test_eft_cancellation_case(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches reassociation or non-nearest behavior that breaks exact cancellation invariants.
    REAL one = A::one();
    REAL neg_one = -one;
    REAL s = one;
    REAL e = one;
    vmplapack::Rtwosum(one, neg_one, s, e);
    require_zero_value(s, tier, "Rtwosum cancellation high part is not zero");
    require_zero_value(e, tier, "Rtwosum cancellation residual is not zero");

    vmplapack::Rfasttwosum(one, neg_one, s, e);
    require_zero_value(s, tier, "Rfasttwosum cancellation high part is not zero");
    require_zero_value(e, tier, "Rfasttwosum cancellation residual is not zero");

    REAL p = A::zero();
    REAL r = one;
    vmplapack::Rtwoproduct(one, neg_one, p, r);
    require<REAL>(p == neg_one, tier, "Rtwoproduct cancellation-sign product is wrong");
    require_zero_value(r, tier, "Rtwoproduct cancellation-sign residual is not zero");
}

template <class REAL>
void test_eft_large_small_case(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches unsafe reassociation that drops the residual in a large-plus-small tie.
    REAL large = tier_power_of_two<REAL>(A::precision_bits());
    REAL s = A::zero();
    REAL e = A::zero();
    vmplapack::Rtwosum(large, A::one(), s, e);
    require<REAL>(s == large, tier, "Rtwosum large/small high part is wrong");
    require<REAL>(e == A::one(), tier, "Rtwosum large/small residual is wrong");

    vmplapack::Rfasttwosum(large, A::one(), s, e);
    require<REAL>(s == large, tier, "Rfasttwosum large/small high part is wrong");
    require<REAL>(e == A::one(), tier, "Rfasttwosum large/small residual is wrong");

    // Catches product-residual breakage on an exact power-of-two product far from unit scale.
    REAL p = A::zero();
    REAL r = A::one();
    vmplapack::Rtwoproduct(large, A::unit_roundoff(), p, r);
    require<REAL>(p == A::one(), tier, "Rtwoproduct large/small product is wrong");
    require_zero_value(r, tier, "Rtwoproduct large/small residual is not zero");
}

template <class REAL>
void test_eft_signed_zero_case(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches implementations that erase the signed-zero product required by IEEE arithmetic.
    REAL zero = A::zero();
    REAL negative_zero = -zero;
    REAL s = A::one();
    REAL e = A::one();
    vmplapack::Rtwosum(negative_zero, negative_zero, s, e);
    require_zero_value(s, tier, "Rtwosum signed-zero high part is not zero");
    require_zero_value(e, tier, "Rtwosum signed-zero residual is not zero");
    require<REAL>(is_negative_zero(s), tier, "Rtwosum did not preserve negative signed zero high part");

    REAL p = A::one();
    REAL r = A::one();
    vmplapack::Rtwoproduct(negative_zero, A::one(), p, r);
    require_zero_value(p, tier, "Rtwoproduct signed-zero product is not zero");
    require_zero_value(r, tier, "Rtwoproduct signed-zero residual is not zero");
    require<REAL>(is_negative_zero(p), tier, "Rtwoproduct did not preserve negative signed zero product");
}

template <class REAL>
void test_eft_tier_contract(const char* tier) {
    test_eft_normal_case<REAL>(tier);
    test_eft_cancellation_case<REAL>(tier);
    test_eft_large_small_case<REAL>(tier);
    test_eft_signed_zero_case<REAL>(tier);
}

template <class REAL>
void test_native_subnormal_eft_case(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches FTZ/DAZ and EFT breakage on additions entirely in the subnormal range.
    REAL eta = A::eta();
    REAL eta2 = eta + eta;
    REAL s = A::zero();
    REAL e = A::one();
    vmplapack::Rtwosum(eta, eta, s, e);
    require<REAL>(s == eta2, tier, "Rtwosum subnormal high part is wrong");
    require_zero_value(e, tier, "Rtwosum subnormal residual is not zero");

    vmplapack::Rfasttwosum(eta, eta, s, e);
    require<REAL>(s == eta2, tier, "Rfasttwosum subnormal high part is wrong");
    require_zero_value(e, tier, "Rfasttwosum subnormal residual is not zero");

    // TwoProduct exactness is tested on normal inputs; Dekker splitting is not required to make
    // exact residual claims for subnormal operands beyond the SPEC §6.2 preconditions.
}

template <class REAL>
void test_native_subnormals_alive(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches FTZ/DAZ settings through an ordinary operation, independent of MXCSR inspection.
    require<REAL>(std::fesetround(FE_TONEAREST) == 0, tier, "failed to set FE_TONEAREST");
    REAL smallest_normal = std::numeric_limits<REAL>::min();
    REAL subnormal = runtime_mul(smallest_normal, A::half());
    require<REAL>(subnormal > A::zero(), tier, "subnormal arithmetic was flushed to zero");
    require<REAL>(subnormal < smallest_normal, tier, "subnormal arithmetic did not produce a subnormal");
}

void test_x86_ftz_daz_bits() {
#if (defined(__i386__) || defined(__x86_64__)) && defined(__SSE__)
    // Catches x86 MXCSR modes that flush subnormal outputs or inputs even when compiler flags look strict.
    unsigned int csr = _mm_getcsr();
    require<double>((csr & (1U << 15U)) == 0U, "x86", "MXCSR FTZ bit is enabled");
    require<double>((csr & (1U << 6U)) == 0U, "x86", "MXCSR DAZ bit is enabled");
#endif
}

void test_native_directed_rounding_runtime() {
    // Catches constant-folding or ignored fenv: upward must round 1 + 2^-53 above 1 at runtime.
    require<double>(std::fesetround(FE_TONEAREST) == 0, "double", "failed to set FE_TONEAREST");
    double half_ulp = 0x1p-53;
    double nearest = runtime_add(1.0, half_ulp);
    require<double>(nearest == 1.0, "double", "nearest runtime addition did not tie to 1");

    double down = 0.0;
    double up = 0.0;
    {
        vmplapack::Rarith<double>::round_down scope;
        down = runtime_add(1.0, half_ulp);
    }
    {
        vmplapack::Rarith<double>::round_up scope;
        up = runtime_add(1.0, half_ulp);
    }
    require<double>(down == 1.0, "double", "downward runtime addition did not stay at 1");
    require<double>(up > 1.0, "double", "upward runtime addition did not rise above 1");

    // Catches multiplication rounding-mode ignorance and accidental contraction-like rewrites.
    double next = std::nextafter(1.0, std::numeric_limits<double>::infinity());
    double mul_down = 0.0;
    double mul_up = 0.0;
    {
        vmplapack::Rarith<double>::round_down scope;
        mul_down = runtime_mul(next, next);
    }
    {
        vmplapack::Rarith<double>::round_up scope;
        mul_up = runtime_mul(next, next);
    }
    require<double>(mul_down < mul_up, "double", "directed runtime multiplication did not change");
    require<double>(std::fegetround() == FE_TONEAREST, "double", "native rounding mode was not restored");
}

void test_native_contract() {
    test_x86_ftz_daz_bits();
    test_native_directed_rounding_runtime();
    test_native_subnormals_alive<float>("float");
    test_native_subnormals_alive<double>("double");
    test_eft_tier_contract<float>("float");
    test_eft_tier_contract<double>("double");
    test_native_subnormal_eft_case<float>("float");
    test_native_subnormal_eft_case<double>("double");
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

void test_mpfr_directed_rounding_runtime(long precision) {
    using REAL = mpfrxx::mpfr_class;
    using A = vmplapack::Rarith<REAL>;

    mpfrxx::set_default_precision_bits(static_cast<mpfr_prec_t>(precision));
    mpfr_set_default_rounding_mode(MPFR_RNDN);

    // Catches an MPFR rounding scope that does not affect ordinary expression evaluation.
    REAL nearest = A::one() + A::unit_roundoff();
    REAL down = A::zero();
    REAL up = A::zero();
    {
        typename A::round_down scope;
        down = A::one() + A::unit_roundoff();
    }
    {
        typename A::round_up scope;
        up = A::one() + A::unit_roundoff();
    }
    require<REAL>(nearest == A::one(), "mpfr", "nearest MPFR addition did not tie to 1");
    require<REAL>(down == A::one(), "mpfr", "downward MPFR addition did not stay at 1");
    require<REAL>(up > A::one(), "mpfr", "upward MPFR addition did not rise above 1");

    // Catches multiplication rounding-mode ignorance in the MPFR tier.
    REAL next = next_after_one_by_unit_roundoff<REAL>();
    REAL mul_down = A::zero();
    REAL mul_up = A::zero();
    {
        typename A::round_down scope;
        mul_down = next * next;
    }
    {
        typename A::round_up scope;
        mul_up = next * next;
    }
    require<REAL>(mul_down < mul_up, "mpfr", "directed MPFR multiplication did not change");
    mpfr_set_default_rounding_mode(MPFR_RNDN);
}

void test_mpfr_contract(long precision) {
    if (precision == 0) {
        return;
    }
    test_mpfr_directed_rounding_runtime(precision);
    test_eft_tier_contract<mpfrxx::mpfr_class>("mpfr");
}
#endif

} // namespace

#undef VMPLAPACK_M9_NOINLINE

int main() {
    test_native_contract();
#ifdef VMPLAPACK_ENABLE_MPFR
    long precision = parse_expected_mpfr_precision();
    test_mpfr_contract(precision);
#endif
    return EXIT_SUCCESS;
}
