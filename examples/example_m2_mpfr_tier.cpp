// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>

namespace {

template <class REAL>
int output_digits() {
    return std::numeric_limits<REAL>::max_digits10;
}

template <>
int output_digits<mpfrxx::mpfr_class>() {
    return 180;
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
void show_tier(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    std::cout << std::setprecision(output_digits<REAL>());

    REAL one = A::one();
    REAL u = A::unit_roundoff();
    REAL a = next_after_one<REAL>();
    REAL b = a;
    REAL p = A::zero();
    REAL e = A::zero();
    vmplapack::Rtwoproduct(a, b, p, e);

    REAL neg_p = -p;
    REAL fma_residual = A::fma(a, b, neg_p);
    REAL recombined = p + e;

    std::cout << tier << '\n';
    std::cout << "  precision_bits = " << A::precision_bits() << '\n';
    std::cout << "  unit_roundoff = " << u << '\n';
    std::cout << "  one = " << one << '\n';
    std::cout << "  next_after_one = " << a << '\n';
    std::cout << "  twoproduct high = " << p << '\n';
    std::cout << "  twoproduct residual = " << e << '\n';
    std::cout << "  fma residual = " << fma_residual << '\n';
    std::cout << "  rounded high + residual = " << recombined << '\n';
}

} // namespace

int main() {
    show_tier<float>("float");
    show_tier<double>("double");
    mpfrxx::set_default_precision_bits(512);
    show_tier<mpfrxx::mpfr_class>("mpfr@512");
    return 0;
}
