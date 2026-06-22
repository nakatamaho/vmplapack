// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>

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
#define VMPLAPACK_TEST_NOINLINE __attribute__((noinline))
#else
#define VMPLAPACK_TEST_NOINLINE
#endif

template <class REAL>
void require(bool condition, const char* tier, const char* message) {
    if (!condition) {
        std::cerr << tier << ": " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

template <class REAL>
VMPLAPACK_TEST_NOINLINE REAL runtime_add(REAL a, REAL b) {
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
VMPLAPACK_TEST_NOINLINE REAL runtime_mul(REAL a, REAL b) {
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
void local_twosum(REAL a, REAL b, REAL& s, REAL& e) {
    s = a + b;
    REAL z = s - a;
    REAL t1 = s - z;
    REAL t2 = a - t1;
    REAL t3 = b - z;
    e = t2 + t3;
}

template <class REAL>
void local_twoproduct(REAL a, REAL b, REAL& p, REAL& e) {
    p = a * b;
    REAL neg_p = -p;
    e = vmplapack::Rarith<REAL>::fma(a, b, neg_p);
}

template <class REAL>
REAL directed_sum_up(REAL a, REAL b, REAL c) {
    using A = vmplapack::Rarith<REAL>;
    typename A::round_up g;
    REAL acc = A::zero();
    REAL t0 = runtime_add(acc, a);
    acc = t0;
    REAL t1 = runtime_add(acc, b);
    acc = t1;
    REAL t2 = runtime_add(acc, c);
    acc = t2;
    return acc;
}

template <class REAL>
REAL directed_sum_down(REAL a, REAL b, REAL c) {
    using A = vmplapack::Rarith<REAL>;
    typename A::round_down g;
    REAL acc = A::zero();
    REAL t0 = runtime_add(acc, a);
    acc = t0;
    REAL t1 = runtime_add(acc, b);
    acc = t1;
    REAL t2 = runtime_add(acc, c);
    acc = t2;
    return acc;
}

template <class REAL>
void test_scope_restoration(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches a copied or leaking upward scope that would leave the process in a directed mode.
    require<REAL>(std::fesetround(FE_TONEAREST) == 0, tier, "failed to set FE_TONEAREST");
    {
        typename A::round_up g;
        require<REAL>(std::fegetround() == FE_UPWARD, tier, "round_up did not set FE_UPWARD");
    }
    require<REAL>(std::fegetround() == FE_TONEAREST, tier, "round_up did not restore FE_TONEAREST");

    // Catches exception-unwind bugs in the downward scope destructor.
    try {
        typename A::round_down g;
        require<REAL>(std::fegetround() == FE_DOWNWARD, tier, "round_down did not set FE_DOWNWARD");
        throw std::runtime_error("scope unwind");
    } catch (const std::runtime_error&) {
    }
    require<REAL>(std::fegetround() == FE_TONEAREST, tier, "round_down did not restore after exception");
}

template <class REAL>
void test_crafted_eft_reconstruction(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches loss of round-to-nearest assumptions in TwoSum on a half-ulp tie at 1.
    require<REAL>(std::fesetround(FE_TONEAREST) == 0, tier, "failed to set FE_TONEAREST");
    REAL s = A::zero();
    REAL e = A::zero();
    local_twosum(A::one(), A::unit_roundoff(), s, e);
    require<REAL>(s == A::one(), tier, "TwoSum crafted high part is wrong");
    require<REAL>(e == A::unit_roundoff(), tier, "TwoSum crafted low part is wrong");

    // Catches a non-correctly-rounded FMA residual in a product with a known power-of-two error.
    REAL next = std::nextafter(A::one(), A::infinity());
    REAL eps = next - A::one();
    REAL p = A::zero();
    REAL r = A::zero();
    local_twoproduct(next, next, p, r);
    REAL eps2 = eps + eps;
    REAL expected_p = A::one() + eps2;
    REAL expected_r = eps * eps;
    require<REAL>(p == expected_p, tier, "TwoProduct crafted high part is wrong");
    require<REAL>(r == expected_r, tier, "TwoProduct crafted low part is wrong");
}

template <class REAL>
void test_subnormals_and_runtime_rounding(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches FTZ/DAZ settings by producing a subnormal through an ordinary operation.
    require<REAL>(std::fesetround(FE_TONEAREST) == 0, tier, "failed to set FE_TONEAREST");
    REAL smallest_normal = std::numeric_limits<REAL>::min();
    REAL subnormal = runtime_mul(smallest_normal, A::half());
    require<REAL>(subnormal > A::zero(), tier, "subnormal result was flushed to zero");
    require<REAL>(subnormal < smallest_normal, tier, "subnormal result was not below min normal");

    // Catches optimization that ignores runtime rounding for a normal addition.
    REAL add_down = A::zero();
    REAL add_up = A::zero();
    {
        typename A::round_down g;
        add_down = runtime_add(A::one(), A::unit_roundoff());
    }
    {
        typename A::round_up g;
        add_up = runtime_add(A::one(), A::unit_roundoff());
    }
    require<REAL>(add_down < add_up, tier, "directed rounding did not change runtime addition");

    // Catches optimization or contraction behavior that ignores runtime rounding for multiplication.
    REAL next = std::nextafter(A::one(), A::infinity());
    REAL mul_down = A::zero();
    REAL mul_up = A::zero();
    {
        typename A::round_down g;
        mul_down = runtime_mul(next, next);
    }
    {
        typename A::round_up g;
        mul_up = runtime_mul(next, next);
    }
    require<REAL>(mul_down < mul_up, tier, "directed rounding did not change runtime multiplication");
    require<REAL>(std::fegetround() == FE_TONEAREST, tier, "rounding mode was not restored");
}

template <class REAL>
void test_directed_pass_brackets_hard_case(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches directed passes that do not enclose a cancellation case that loses 1 in nearest mode.
    require<REAL>(std::fesetround(FE_TONEAREST) == 0, tier, "failed to set FE_TONEAREST");
    REAL large = std::ldexp(A::one(), static_cast<int>(A::precision_bits()));
    REAL lo = directed_sum_down(large, A::one(), -large);
    REAL hi = directed_sum_up(large, A::one(), -large);
    require<REAL>(lo <= A::one(), tier, "directed lower pass is above the true value");
    require<REAL>(A::one() <= hi, tier, "directed upper pass is below the true value");
    require<REAL>(lo <= hi, tier, "directed passes produced inverted bounds");
}

template <class REAL>
void test_tier(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    static_assert(std::is_same<decltype(A::precision_bits()), long>::value,
                  "precision_bits must return signed long");
    static_assert(!std::is_copy_constructible<typename A::round_up>::value,
                  "round_up must be non-copyable");
    static_assert(!std::is_move_constructible<typename A::round_up>::value,
                  "round_up must be non-movable");
    static_assert(!std::is_copy_constructible<typename A::round_down>::value,
                  "round_down must be non-copyable");
    static_assert(!std::is_move_constructible<typename A::round_down>::value,
                  "round_down must be non-movable");

    test_scope_restoration<REAL>(tier);
    test_crafted_eft_reconstruction<REAL>(tier);
    test_subnormals_and_runtime_rounding<REAL>(tier);
    test_directed_pass_brackets_hard_case<REAL>(tier);
}

} // namespace

#undef VMPLAPACK_TEST_NOINLINE

int main() {
    test_tier<float>("float");
    test_tier<double>("double");
    return EXIT_SUCCESS;
}
