// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include "../tests/Rgendot.h"

#include <vmplapack/vmplapack.h>

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>

namespace {

template <class REAL>
int output_digits() {
    return std::numeric_limits<REAL>::max_digits10;
}

template <>
int output_digits<mpfrxx::mpfr_class>() {
    return 90;
}

const char* status_name(vmplapack::Rstatus status) {
    switch (status) {
        case vmplapack::Rstatus::ok:
            return "ok";
        case vmplapack::Rstatus::unbounded:
            return "unbounded";
        case vmplapack::Rstatus::non_finite:
            return "non_finite";
        case vmplapack::Rstatus::invalid_input:
            return "invalid_input";
    }
    return "invalid_input";
}

template <class REAL>
REAL naive_dot(const vmplapack::gendot::Rdot_case<REAL>& c) {
    using A = vmplapack::Rarith<REAL>;

    REAL acc = A::zero();
    for (std::size_t i = 0; i < c.x.size(); ++i) {
        REAL prod = c.x[i] * c.y[i];
        REAL next = acc + prod;
        acc = next;
    }
    return acc;
}

template <class REAL>
REAL directed_dot_up(const vmplapack::gendot::Rdot_case<REAL>& c) {
    using A = vmplapack::Rarith<REAL>;

    typename A::round_up scope;
    REAL acc = A::zero();
    for (std::size_t i = 0; i < c.x.size(); ++i) {
        REAL prod = c.x[i] * c.y[i];
        REAL next = acc + prod;
        acc = next;
    }
    return acc;
}

template <class REAL>
REAL directed_dot_down(const vmplapack::gendot::Rdot_case<REAL>& c) {
    using A = vmplapack::Rarith<REAL>;

    typename A::round_down scope;
    REAL acc = A::zero();
    for (std::size_t i = 0; i < c.x.size(); ++i) {
        REAL prod = c.x[i] * c.y[i];
        REAL next = acc + prod;
        acc = next;
    }
    return acc;
}

template <class REAL>
void show_terms(const vmplapack::gendot::Rdot_case<REAL>& c) {
    using A = vmplapack::Rarith<REAL>;

    REAL prefix = A::zero();
    std::cout << "  terms under ordinary nearest rounding" << '\n';
    for (std::size_t i = 0; i < c.x.size(); ++i) {
        REAL prod = c.x[i] * c.y[i];
        REAL next = prefix + prod;
        prefix = next;
        std::cout << "    i=" << i
                  << " x=" << c.x[i]
                  << " y=" << c.y[i]
                  << " x*y=" << prod
                  << " prefix=" << prefix << '\n';
    }
}

template <class REAL>
void show_box(const vmplapack::Rmidrad<REAL>& box) {
    REAL lo = box.mid - box.rad;
    REAL hi = box.mid + box.rad;
    std::cout << "  verified mid = " << box.mid << '\n';
    std::cout << "  verified rad = " << box.rad << '\n';
    std::cout << "  displayed enclosure lower = " << lo << '\n';
    std::cout << "  displayed enclosure upper = " << hi << '\n';
    std::cout << "  status = " << status_name(box.status) << '\n';
    std::cout << "  status_rank = " << vmplapack::Rstatus_rank(box.status) << '\n';
}

template <class REAL>
void show_tier(const char* tier, int exponent) {
    using A = vmplapack::Rarith<REAL>;

    vmplapack::gendot::Rdot_case<REAL> c =
        vmplapack::gendot::family_c_exponent_cancellation<REAL>(exponent, A::one());
    std::ptrdiff_t n = static_cast<std::ptrdiff_t>(c.x.size());
    REAL naive = naive_dot(c);
    REAL accurate = vmplapack::Rdot(n, c.x.data(), 1, c.y.data(), 1);
    REAL lower_pass = directed_dot_down(c);
    REAL upper_pass = directed_dot_up(c);
    vmplapack::Rmidrad<REAL> box = vmplapack::vRdot(n, c.x.data(), 1, c.y.data(), 1);

    std::cout << std::setprecision(output_digits<REAL>());
    std::cout << tier << '\n';
    std::cout << "  case = [2^" << exponent << ", 1, 2^" << exponent << "] dot [1, 1, -1]" << '\n';
    std::cout << "  exact = " << c.exact << '\n';
    show_terms(c);
    std::cout << "  naive dot = " << naive << '\n';
    std::cout << "  accurate Rdot = " << accurate << '\n';
    std::cout << "  directed lower pass = " << lower_pass << '\n';
    std::cout << "  directed upper pass = " << upper_pass << '\n';
    show_box(box);
    std::cout << "  exact covered by displayed enclosure = "
              << std::boolalpha << ((box.mid - box.rad) <= c.exact && c.exact <= (box.mid + box.rad)) << '\n';
}

} // namespace

int main() {
    show_tier<float>("float", 30);
    show_tier<double>("double", 60);
    mpfrxx::set_default_precision_bits(512);
    show_tier<mpfrxx::mpfr_class>("mpfr@512", 520);
    return 0;
}
