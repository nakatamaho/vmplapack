// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include "../tests/Rdot_oracle.h"
#include "../tests/Rgendot.h"

#include <vmplapack/vmplapack.h>

#include <cstddef>
#include <cstdlib>
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

mpfrxx::mpfr_class upward_abs_diff(const mpfrxx::mpfr_class& a,
                                   const mpfrxx::mpfr_class& b,
                                   mpfr_prec_t precision) {
    mpfrxx::mpfr_class aa = vmplapack::oracle::widen_mpfr(a, precision);
    mpfrxx::mpfr_class bb = vmplapack::oracle::widen_mpfr(b, precision);
    mpfrxx::mpfr_class result = mpfrxx::mpfr_class::with_precision(precision);
    mpfrxx::rounding_mode_scope scope(MPFR_RNDU);
    if (aa >= bb) {
        result = aa - bb;
    } else {
        result = bb - aa;
    }
    return result;
}

mpfrxx::mpfr_class relative_from_abs(const mpfrxx::mpfr_class& abs_value,
                                     const mpfrxx::mpfr_class& oracle_mid,
                                     mpfr_prec_t precision) {
    mpfrxx::mpfr_class abs_mid = vmplapack::oracle::abs_at(oracle_mid, precision);
    if (abs_mid == vmplapack::oracle::zero_at(precision)) {
        if (abs_value == vmplapack::oracle::zero_at(precision)) {
            return vmplapack::oracle::zero_at(precision);
        }
        return vmplapack::gendot::infinity_at(precision);
    }
    mpfrxx::rounding_mode_scope scope(MPFR_RNDU);
    mpfrxx::mpfr_class relative = abs_value / abs_mid;
    return relative;
}

template <class REAL>
mpfrxx::mpfr_class absolute_error(const REAL& value,
                                  const mpfrxx::mpfr_class& oracle_mid,
                                  mpfr_prec_t precision) {
    mpfrxx::mpfr_class value_mp = vmplapack::oracle::widen_value(value, precision);
    return upward_abs_diff(value_mp, oracle_mid, precision);
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
void print_scalar_result(const char* label,
                         const REAL& value,
                         const mpfrxx::mpfr_class& oracle_mid,
                         mpfr_prec_t precision) {
    mpfrxx::mpfr_class err = absolute_error(value, oracle_mid, precision);
    mpfrxx::mpfr_class rel = relative_from_abs(err, oracle_mid, precision);
    std::cout << "  " << label << " = " << value << '\n';
    std::cout << "    abs error = " << err << '\n';
    std::cout << "    rel error = " << rel << '\n';
    std::cout << "    radius = null" << '\n';
    std::cout << "    enclosed = null" << '\n';
}

template <class REAL>
void print_box_result(const char* label,
                      const vmplapack::Rmidrad<REAL>& box,
                      const vmplapack::oracle::Rdot_interval& ref,
                      const mpfrxx::mpfr_class& oracle_mid,
                      mpfr_prec_t precision) {
    mpfrxx::mpfr_class err = absolute_error(box.mid, oracle_mid, precision);
    mpfrxx::mpfr_class rel = relative_from_abs(err, oracle_mid, precision);
    mpfrxx::mpfr_class rad = vmplapack::oracle::widen_value(box.rad, precision);
    mpfrxx::mpfr_class rel_rad = relative_from_abs(rad, oracle_mid, precision);
    std::cout << "  " << label << '\n';
    std::cout << "    mid = " << box.mid << '\n';
    std::cout << "    rad = " << box.rad << '\n';
    std::cout << "    relative radius = " << rel_rad << '\n';
    std::cout << "    mid abs error = " << err << '\n';
    std::cout << "    mid rel error = " << rel << '\n';
    std::cout << "    status = " << status_name(box.status) << '\n';
    std::cout << "    enclosed oracle interval = " << std::boolalpha << midrad_covers_oracle(box, ref) << '\n';
}

template <class REAL>
void show_tier(const char* tier, int target_power2) {
    vmplapack::gendot::Rgenerated_dot_case<REAL> generated =
        vmplapack::gendot::randomized_high_condition_power2<REAL>(target_power2,
                                                                  11U,
                                                                  vmplapack::gendot::Rpermutation::sorted);
    if (generated.status != vmplapack::gendot::Rgenerator_status::ok) {
        std::cerr << tier << ": generator did not produce an ok case\n";
        std::exit(EXIT_FAILURE);
    }

    const vmplapack::gendot::Rdot_case<REAL>& c = generated.data;
    std::ptrdiff_t n = static_cast<std::ptrdiff_t>(c.x.size());
    vmplapack::oracle::Rdot_interval ref = vmplapack::oracle::Rdot_oracle(n, c.x.data(), 1, c.y.data(), 1);
    mpfr_prec_t precision = ref.precision + static_cast<mpfr_prec_t>(vmplapack::Rarith<REAL>::precision_bits() + 192L);
    mpfrxx::mpfr_class oracle_mid = interval_midpoint(ref, precision);

    REAL naive = naive_dot(c);
    REAL accurate = vmplapack::Rdot(n, c.x.data(), 1, c.y.data(), 1);
    vmplapack::Rmidrad<REAL> directed = vmplapack::vRdot(n, c.x.data(), 1, c.y.data(), 1);
    vmplapack::Rmidrad<REAL> apriori = vmplapack::vRdot_apriori(n, c.x.data(), 1, c.y.data(), 1);

    std::cout << std::setprecision(output_digits<REAL>()) << std::boolalpha;
    std::cout << tier << '\n';
    std::cout << "  generator = m8_random_power2" << '\n';
    std::cout << "  seed = " << generated.seed << '\n';
    std::cout << "  target cond_oro = " << generated.target_cond_oro << '\n';
    std::cout << "  realized cond_oro = " << generated.measured_cond_oro << '\n';
    std::cout << "  realized cond_sum = " << generated.measured_cond_sum << '\n';
    std::cout << "  n = " << n << '\n';
    std::cout << "  oracle lo = " << ref.lo << '\n';
    std::cout << "  oracle hi = " << ref.hi << '\n';
    std::cout << "  oracle midpoint = " << oracle_mid << '\n';
    print_scalar_result("naive dot", naive, oracle_mid, precision);
    print_scalar_result("Rdot", accurate, oracle_mid, precision);
    print_box_result("vRdot", directed, ref, oracle_mid, precision);
    print_box_result("vRdot_apriori", apriori, ref, oracle_mid, precision);
}

} // namespace

int main() {
    show_tier<float>("float", 20);
    show_tier<double>("double", 40);
    mpfrxx::set_default_precision_bits(512);
    show_tier<mpfrxx::mpfr_class>("mpfr@512", 128);
    return EXIT_SUCCESS;
}
