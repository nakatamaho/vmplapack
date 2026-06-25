// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

#include <cstddef>
#include <cstdlib>
#include <iostream>

namespace {

template <class REAL>
int check_box(const vmplapack::Rmidrad<REAL>& box, REAL expected, const char* tier) {
    if (box.status != vmplapack::Rstatus::ok) {
        std::cerr << tier << ": verified dot did not return ok\n";
        return EXIT_FAILURE;
    }
    if (box.rad < vmplapack::Rarith<REAL>::zero()) {
        std::cerr << tier << ": verified dot returned a negative radius\n";
        return EXIT_FAILURE;
    }
    REAL lo = box.mid - box.rad;
    REAL hi = box.mid + box.rad;
    if (!(lo <= expected && expected <= hi)) {
        std::cerr << tier << ": verified dot does not enclose the expected value\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int check_double_dot() {
    double x[3] = {1.0, 2.0, 3.0};
    double y[3] = {4.0, 5.0, 6.0};
    vmplapack::Rmidrad<double> box = vmplapack::vRdot(3, x, 1, y, 1);
    return check_box(box, 32.0, "double");
}

#ifdef VMPLAPACK_CONSUMER_EXPECT_MPFR
int check_mpfr_dot() {
#ifndef VMPLAPACK_ENABLE_MPFR
#error "The consumer expected MPFR support, but the installed vMPLAPACK target did not define VMPLAPACK_ENABLE_MPFR."
#endif
    mpfrxx::set_default_precision_bits(512);
    mpfrxx::mpfr_class x[3] = {
        mpfrxx::mpfr_class::with_precision(512, 1.0),
        mpfrxx::mpfr_class::with_precision(512, 2.0),
        mpfrxx::mpfr_class::with_precision(512, 3.0)};
    mpfrxx::mpfr_class y[3] = {
        mpfrxx::mpfr_class::with_precision(512, 4.0),
        mpfrxx::mpfr_class::with_precision(512, 5.0),
        mpfrxx::mpfr_class::with_precision(512, 6.0)};
    mpfrxx::mpfr_class expected = mpfrxx::mpfr_class::with_precision(512, 32.0);
    vmplapack::Rmidrad<mpfrxx::mpfr_class> box = vmplapack::vRdot(3, x, 1, y, 1);
    return check_box(box, expected, "mpfr");
}
#endif

} // namespace

int main() {
    int double_status = check_double_dot();
    if (double_status != EXIT_SUCCESS) {
        return double_status;
    }
#ifdef VMPLAPACK_CONSUMER_EXPECT_MPFR
    int mpfr_status = check_mpfr_dot();
    if (mpfr_status != EXIT_SUCCESS) {
        return mpfr_status;
    }
#endif
    return EXIT_SUCCESS;
}
