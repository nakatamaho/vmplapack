// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

#include <cfenv>
#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>

namespace {

template <class REAL>
int output_digits() {
    return std::numeric_limits<REAL>::max_digits10;
}

#ifdef VMPLAPACK_ENABLE_MPFR
template <>
int output_digits<mpfrxx::mpfr_class>() {
    return 180;
}
#endif

template <class REAL>
REAL runtime_add(REAL a, REAL b) {
#ifdef VMPLAPACK_ENABLE_MPFR
    if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        REAL z = a + b;
        return z;
    } else
#endif
    {
        volatile REAL va = a;
        volatile REAL vb = b;
        REAL x = va;
        REAL y = vb;
        REAL z = x + y;
        return z;
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
void show_midrad(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    std::cout << std::setprecision(output_digits<REAL>());

    REAL lo = A::one();
    REAL hi = next_after_one<REAL>();
    REAL mid = vmplapack::Rmidpoint(lo, hi);
    vmplapack::Rmidrad<REAL> box = vmplapack::Rmake_midrad(lo, hi, mid);

    std::cout << tier << '\n';
    std::cout << "  precision_bits = " << A::precision_bits() << '\n';
    std::cout << "  lo = " << lo << '\n';
    std::cout << "  hi = " << hi << '\n';
    std::cout << "  mid = " << box.mid << '\n';
    std::cout << "  rad = " << box.rad << '\n';
    std::cout << "  status_rank = " << vmplapack::Rstatus_rank(box.status) << '\n';

    REAL rounded_up = A::zero();
    {
        typename A::round_up scope;
        rounded_up = runtime_add(lo, A::unit_roundoff());
    }

    std::cout << "  upward(1 + u) = " << rounded_up << '\n';
    std::cout << "  rounding_restored = " << (std::fegetround() == FE_TONEAREST) << '\n';
}

} // namespace

int main() {
    show_midrad<float>("float");
    show_midrad<double>("double");
#ifdef VMPLAPACK_ENABLE_MPFR
    mpfrxx::set_default_precision_bits(512);
    show_midrad<mpfrxx::mpfr_class>("mpfr@512");
#endif
    return 0;
}
