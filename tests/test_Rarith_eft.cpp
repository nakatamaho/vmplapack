// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

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

template <class REAL>
void require(bool condition, const char* tier, const char* message) {
    if (!condition) {
        std::cerr << tier << ": " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

template <class REAL>
REAL next_after_one() {
    using A = vmplapack::Rarith<REAL>;

    REAL two = A::one() + A::one();
    REAL step = two * A::unit_roundoff();
    REAL next = A::one() + step;
    return next;
}

template <class REAL>
void test_twosum_crafted(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches a TwoSum implementation that is not exact on a nearest-even half-ulp case.
    require<REAL>(std::fesetround(FE_TONEAREST) == 0, tier, "failed to set FE_TONEAREST");
    REAL s = A::zero();
    REAL e = A::zero();
    vmplapack::Rtwosum(A::one(), A::unit_roundoff(), s, e);
    require<REAL>(s == A::one(), tier, "Rtwosum high part is wrong");
    require<REAL>(e == A::unit_roundoff(), tier, "Rtwosum low part is wrong");

    // Catches input/output aliasing bugs in calls shaped like Rtwosum(s, x, s, e).
    s = A::one();
    e = A::zero();
    vmplapack::Rtwosum(s, A::unit_roundoff(), s, e);
    require<REAL>(s == A::one(), tier, "aliased Rtwosum high part is wrong");
    require<REAL>(e == A::unit_roundoff(), tier, "aliased Rtwosum low part is wrong");

    // Catches input/output aliasing bugs in calls shaped like Rtwosum(p, h, p, q).
    REAL p = A::one();
    REAL q = A::zero();
    vmplapack::Rtwosum(p, A::unit_roundoff(), p, q);
    require<REAL>(p == A::one(), tier, "dot-shaped aliased Rtwosum high part is wrong");
    require<REAL>(q == A::unit_roundoff(), tier, "dot-shaped aliased Rtwosum low part is wrong");
}

template <class REAL>
void test_fasttwosum_crafted(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches a FastTwoSum implementation that violates the abs(a) >= abs(b) exact case.
    require<REAL>(std::fesetround(FE_TONEAREST) == 0, tier, "failed to set FE_TONEAREST");
    REAL s = A::zero();
    REAL e = A::zero();
    vmplapack::Rfasttwosum(A::one(), A::unit_roundoff(), s, e);
    require<REAL>(s == A::one(), tier, "Rfasttwosum high part is wrong");
    require<REAL>(e == A::unit_roundoff(), tier, "Rfasttwosum low part is wrong");

    // Catches input/output aliasing bugs in FastTwoSum.
    s = A::one();
    e = A::zero();
    vmplapack::Rfasttwosum(s, A::unit_roundoff(), s, e);
    require<REAL>(s == A::one(), tier, "aliased Rfasttwosum high part is wrong");
    require<REAL>(e == A::unit_roundoff(), tier, "aliased Rfasttwosum low part is wrong");
}

template <class REAL>
void test_split_basic(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches a broken split factor for values that should split without a low part.
    REAL hi = A::zero();
    REAL lo = A::one();
    REAL three = A::one() + A::one();
    three = three + A::one();
    REAL value = three * A::half();
    vmplapack::Rsplit(value, hi, lo);
    require<REAL>(hi == value, tier, "Rsplit high part is wrong for an exactly high value");
    require<REAL>(lo == A::zero(), tier, "Rsplit low part is wrong for an exactly high value");
}

template <class REAL>
void test_twoproduct_crafted(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches a TwoProduct implementation whose high or low part disagrees with a known product.
    require<REAL>(std::fesetround(FE_TONEAREST) == 0, tier, "failed to set FE_TONEAREST");
    REAL next = next_after_one<REAL>();
    REAL eps = next - A::one();
    REAL p = A::zero();
    REAL e = A::zero();
    vmplapack::Rtwoproduct(next, next, p, e);
    REAL eps2 = eps + eps;
    REAL expected_p = A::one() + eps2;
    REAL expected_e = eps * eps;
    require<REAL>(p == expected_p, tier, "Rtwoproduct high part is wrong");
    require<REAL>(e == expected_e, tier, "Rtwoproduct low part is wrong");

    // Catches input/output aliasing bugs when the product high part aliases an input.
    REAL aliased = next;
    e = A::zero();
    vmplapack::Rtwoproduct(aliased, next, aliased, e);
    require<REAL>(aliased == expected_p, tier, "aliased Rtwoproduct high part is wrong");
    require<REAL>(e == expected_e, tier, "aliased Rtwoproduct low part is wrong");
}

template <class REAL>
void test_tier(const char* tier) {
    test_twosum_crafted<REAL>(tier);
    test_fasttwosum_crafted<REAL>(tier);
    test_split_basic<REAL>(tier);
    test_twoproduct_crafted<REAL>(tier);
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

void require_mpfr_precision(const mpfrxx::mpfr_class& value, long expected, const char* message) {
    if (mpfr_get_prec(value.mpfr_data()) != static_cast<mpfr_prec_t>(expected)) {
        std::cerr << "mpfr: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

void test_mpfr_precision_discipline(long expected) {
    using REAL = mpfrxx::mpfr_class;
    using A = vmplapack::Rarith<REAL>;

    mpfrxx::set_default_precision_bits(static_cast<mpfr_prec_t>(expected));
    require<REAL>(A::precision_bits() == expected, "mpfr", "default precision does not match W");

    REAL zero = A::zero();
    REAL one = A::one();
    REAL half = A::half();
    REAL u = A::unit_roundoff();
    require_mpfr_precision(zero, expected, "zero precision is not W");
    require_mpfr_precision(one, expected, "one precision is not W");
    require_mpfr_precision(half, expected, "half precision is not W");
    require_mpfr_precision(u, expected, "unit_roundoff precision is not W");
}

void test_mpfr_rounding_scope(long expected) {
    using REAL = mpfrxx::mpfr_class;
    using A = vmplapack::Rarith<REAL>;

    mpfrxx::set_default_precision_bits(static_cast<mpfr_prec_t>(expected));
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
    require<REAL>(down < up, "mpfr", "directed rounding did not change MPFR addition");
}
#endif

} // namespace

int main() {
    test_tier<float>("float");
    test_tier<double>("double");
#ifdef VMPLAPACK_ENABLE_MPFR
    long expected = parse_expected_mpfr_precision();
    if (expected != 0) {
        test_mpfr_precision_discipline(expected);
        test_mpfr_rounding_scope(expected);
        test_tier<mpfrxx::mpfr_class>("mpfr");
    }
#endif
    return EXIT_SUCCESS;
}
