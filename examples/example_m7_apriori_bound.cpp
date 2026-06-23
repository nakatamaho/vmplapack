// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include "../tests/Rdot_oracle.h"
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

mpfrxx::mpfr_class interval_midpoint(const vmplapack::oracle::Rdot_interval& interval,
                                     mpfr_prec_t precision) {
    mpfrxx::mpfr_class lo = vmplapack::oracle::widen_mpfr(interval.lo, precision);
    mpfrxx::mpfr_class hi = vmplapack::oracle::widen_mpfr(interval.hi, precision);
    mpfrxx::mpfr_class sum = mpfrxx::mpfr_class::with_precision(precision);
    sum = lo + hi;
    mpfrxx::mpfr_class half = mpfrxx::mpfr_class::with_precision(precision, 0.5);
    mpfrxx::mpfr_class mid = mpfrxx::mpfr_class::with_precision(precision);
    mid = sum * half;
    return mid;
}

template <class REAL>
bool covers_exact_display(const vmplapack::Rmidrad<REAL>& box, const REAL& exact) {
    REAL lo = box.mid - box.rad;
    REAL hi = box.mid + box.rad;
    return lo <= exact && exact <= hi;
}

template <class REAL>
void show_tier(const char* tier, int exponent) {
    using A = vmplapack::Rarith<REAL>;

    vmplapack::gendot::Rdot_case<REAL> c =
        vmplapack::gendot::family_c_exponent_cancellation<REAL>(exponent, A::one());
    std::ptrdiff_t n = static_cast<std::ptrdiff_t>(c.x.size());
    REAL naive = naive_dot(c);
    REAL accurate = vmplapack::Rdot(n, c.x.data(), 1, c.y.data(), 1);
    REAL s_up = vmplapack::Rsum_abs_dot_upper(n, c.x.data(), 1, c.y.data(), 1);
    vmplapack::Rmidrad<REAL> directed = vmplapack::vRdot(n, c.x.data(), 1, c.y.data(), 1);
    vmplapack::Rmidrad<REAL> apriori = vmplapack::vRdot_apriori(n, c.x.data(), 1, c.y.data(), 1);

    vmplapack::oracle::Rdot_interval interval = vmplapack::oracle::Rdot_oracle(n, c.x.data(), 1, c.y.data(), 1);
    mpfr_prec_t precision = interval.precision + 128;
    mpfrxx::mpfr_class oracle_mid = interval_midpoint(interval, precision);
    mpfrxx::mpfr_class condition = vmplapack::oracle::condition_estimate(n,
                                                                         c.x.data(),
                                                                         1,
                                                                         c.y.data(),
                                                                         1,
                                                                         c.exact,
                                                                         precision);

    std::cout << std::setprecision(output_digits<REAL>()) << std::boolalpha;
    std::cout << tier << '\n';
    std::cout << "  case = [2^" << exponent << ", 1, 2^" << exponent << "] dot [1, 1, -1]" << '\n';
    std::cout << "  exact = " << c.exact << '\n';
    std::cout << "  condition estimate = " << condition << '\n';
    std::cout << "  oracle midpoint = " << oracle_mid << '\n';
    std::cout << "  naive dot = " << naive << '\n';
    std::cout << "  accurate Rdot = " << accurate << '\n';
    std::cout << "  S_up = " << s_up << '\n';
    std::cout << "  vRdot directed rad = " << directed.rad << '\n';
    std::cout << "  vRdot directed status = " << status_name(directed.status) << '\n';
    std::cout << "  vRdot_apriori mid = " << apriori.mid << '\n';
    std::cout << "  vRdot_apriori rad = " << apriori.rad << '\n';
    std::cout << "  vRdot_apriori status = " << status_name(apriori.status) << '\n';
    std::cout << "  apriori tighter than directed = " << (apriori.rad < directed.rad) << '\n';
    std::cout << "  exact covered by apriori display = " << covers_exact_display(apriori, c.exact) << '\n';
}

} // namespace

int main() {
    show_tier<float>("float", 30);
    show_tier<double>("double", 60);
    mpfrxx::set_default_precision_bits(512);
    show_tier<mpfrxx::mpfr_class>("mpfr@512", 700);
    return 0;
}
