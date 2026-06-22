// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

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
void show_box(const vmplapack::Rmidrad<REAL>& box, std::size_t row) {
    REAL lo = box.mid - box.rad;
    REAL hi = box.mid + box.rad;
    std::cout << "  row " << row << '\n';
    std::cout << "    mid = " << box.mid << '\n';
    std::cout << "    rad = " << box.rad << '\n';
    std::cout << "    lower display = " << lo << '\n';
    std::cout << "    upper display = " << hi << '\n';
    std::cout << "    status = " << status_name(box.status) << '\n';
}

template <class REAL>
void show_tier(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    REAL two = A::one() + A::one();
    REAL three = two + A::one();
    REAL four = two + two;
    REAL twenty = four * (four + A::one());

    std::vector<REAL> matrix;
    matrix.push_back(two);
    matrix.push_back(A::one());
    matrix.push_back(-A::one());
    matrix.push_back(four);
    std::vector<REAL> x;
    x.push_back(three);
    x.push_back(four);
    std::vector<REAL> b;
    b.push_back(twenty);
    b.push_back(-A::one());
    std::vector<vmplapack::Rmidrad<REAL>> out(static_cast<std::size_t>(2));

    vmplapack::Rstatus status = vmplapack::vRresidual<REAL>(2, 2, matrix.data(), 2, x.data(), b.data(), out.data());

    std::cout << std::setprecision(output_digits<REAL>());
    std::cout << tier << '\n';
    std::cout << "  A = [[2, 1], [-1, 4]]" << '\n';
    std::cout << "  x = [3, 4]" << '\n';
    std::cout << "  b = [20, -1]" << '\n';
    std::cout << "  residual = b - A*x" << '\n';
    std::cout << "  return status = " << status_name(status) << '\n';
    std::cout << "  return status_rank = " << vmplapack::Rstatus_rank(status) << '\n';
    show_box(out[0], 0);
    show_box(out[1], 1);
}

} // namespace

int main() {
    show_tier<float>("float");
    show_tier<double>("double");
    mpfrxx::set_default_precision_bits(512);
    show_tier<mpfrxx::mpfr_class>("mpfr@512");
    return 0;
}
