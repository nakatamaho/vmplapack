// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

#include <iomanip>
#include <iostream>
#include <limits>

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
REAL next_after_one() {
    using A = vmplapack::Rarith<REAL>;

    REAL two = A::one() + A::one();
    REAL step = two * A::unit_roundoff();
    REAL next = A::one() + step;
    return next;
}

template <class REAL>
void show_twoproduct(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    std::cout << std::setprecision(output_digits<REAL>());

    REAL one = A::one();
    REAL a = next_after_one<REAL>();
    REAL b = a;
    REAL p = A::zero();
    REAL e = A::zero();
    vmplapack::Rtwoproduct(a, b, p, e);

    REAL eps = a - one;
    REAL expected_residual = eps * eps;
    REAL rounded_recombine = p + e;

    std::cout << tier << '\n';
    std::cout << "  a = " << a << '\n';
    std::cout << "  rounded product p = " << p << '\n';
    std::cout << "  residual e = " << e << '\n';
    std::cout << "  expected eps^2 = " << expected_residual << '\n';
    std::cout << "  rounded p + e = " << rounded_recombine << '\n';
}

} // namespace

int main() {
    show_twoproduct<float>("float");
    show_twoproduct<double>("double");
#ifdef VMPLAPACK_ENABLE_MPFR
    mpfrxx::set_default_precision_bits(512);
    show_twoproduct<mpfrxx::mpfr_class>("mpfr@512");
#endif
    return 0;
}
