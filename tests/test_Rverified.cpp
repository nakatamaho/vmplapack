// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include "Rdot_oracle.h"
#include "Rgendot.h"

#include <vmplapack/vmplapack.h>

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

namespace {

template <class REAL>
void require(bool condition, const char* tier, const char* message) {
    if (!condition) {
        std::cerr << tier << ": " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

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

template <class REAL>
REAL exact_epsilon_above_one() {
    using A = vmplapack::Rarith<REAL>;

    REAL two = A::one() + A::one();
    REAL step = two * A::unit_roundoff();
    REAL next = A::one() + step;
    REAL eps = next - A::one();
    return eps;
}

template <class REAL>
std::vector<REAL> ones(std::size_t n) {
    std::vector<REAL> y;
    y.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        y.push_back(vmplapack::Rarith<REAL>::one());
    }
    return y;
}

template <class REAL>
REAL quiet_nan_value() {
    if constexpr (std::is_same<REAL, float>::value || std::is_same<REAL, double>::value) {
        return std::numeric_limits<REAL>::quiet_NaN();
    } else if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        REAL value = REAL::with_precision(static_cast<mpfr_prec_t>(vmplapack::Rarith<REAL>::precision_bits()));
        mpfr_set_nan(value.mpfr_data());
        return value;
    } else {
        static_assert(vmplapack::Ralways_false<REAL>::value, "quiet_nan_value does not support this REAL type.");
    }
}

template <class REAL>
REAL overflow_large_value() {
    if constexpr (std::is_same<REAL, float>::value || std::is_same<REAL, double>::value) {
        return std::numeric_limits<REAL>::max();
    } else if constexpr (std::is_same<REAL, mpfrxx::mpfr_class>::value) {
        REAL value = REAL::with_precision(static_cast<mpfr_prec_t>(vmplapack::Rarith<REAL>::precision_bits()));
        mpfr_exp_t exponent = mpfr_get_emax() - 1;
        mpfr_set_ui_2exp(value.mpfr_data(), 1, exponent, MPFR_RNDN);
        return value;
    } else {
        static_assert(vmplapack::Ralways_false<REAL>::value, "overflow_large_value does not support this REAL type.");
    }
}

template <class REAL>
bool midrad_covers_oracle(const vmplapack::Rmidrad<REAL>& m,
                          const vmplapack::oracle::Rdot_interval& ref) {
    using A = vmplapack::Rarith<REAL>;

    if (m.status == vmplapack::Rstatus::unbounded) {
        return !A::is_finite(m.rad);
    }
    if (m.status != vmplapack::Rstatus::ok || m.rad < A::zero()) {
        return false;
    }

    mpfr_prec_t precision = ref.precision + static_cast<mpfr_prec_t>(A::precision_bits() + 128L);
    mpfrxx::mpfr_class mid = vmplapack::oracle::widen_value(m.mid, precision);
    mpfrxx::mpfr_class rad = vmplapack::oracle::widen_value(m.rad, precision);
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
void require_vrdot_contains_oracle(const vmplapack::gendot::Rdot_case<REAL>& c, const char* tier) {
    // Catches a verified dot interval that is close to the oracle midpoint but does not enclose the oracle interval.
    vmplapack::Rmidrad<REAL> m = vmplapack::vRdot(static_cast<std::ptrdiff_t>(c.x.size()),
                                                 c.x.data(),
                                                 1,
                                                 c.y.data(),
                                                 1);
    vmplapack::oracle::Rdot_interval ref =
        vmplapack::oracle::Rdot_oracle(static_cast<std::ptrdiff_t>(c.x.size()), c.x.data(), 1, c.y.data(), 1);
    require<REAL>(m.status == vmplapack::Rstatus::ok, tier, "vRdot returned non-ok for a finite regular case");
    require<REAL>(midrad_covers_oracle(m, ref), tier, "vRdot did not cover the oracle interval");
}

template <class REAL>
void require_vrsum_contains_oracle(const std::vector<REAL>& x, const char* tier) {
    // Catches a verified sum interval that does not enclose the equivalent dot(x, 1) oracle interval.
    std::vector<REAL> y = ones<REAL>(x.size());
    vmplapack::Rmidrad<REAL> m = vmplapack::vRsum(static_cast<std::ptrdiff_t>(x.size()), x.data(), 1);
    vmplapack::oracle::Rdot_interval ref =
        vmplapack::oracle::Rdot_oracle(static_cast<std::ptrdiff_t>(x.size()), x.data(), 1, y.data(), 1);
    require<REAL>(m.status == vmplapack::Rstatus::ok, tier, "vRsum returned non-ok for a finite regular case");
    require<REAL>(midrad_covers_oracle(m, ref), tier, "vRsum did not cover the oracle interval");
}

template <class REAL>
std::vector<REAL> alternating_sum_case(int pairs, REAL delta) {
    using A = vmplapack::Rarith<REAL>;

    std::vector<REAL> x;
    x.reserve(static_cast<std::size_t>(2 * pairs + 1));
    for (int i = 0; i < pairs; ++i) {
        x.push_back(A::one());
        x.push_back(-A::one());
    }
    x.push_back(delta);
    return x;
}

template <class REAL>
std::vector<REAL> exponent_sum_case(int exponent) {
    using A = vmplapack::Rarith<REAL>;

    std::vector<REAL> x;
    REAL large = vmplapack::gendot::power_of_two<REAL>(exponent);
    x.push_back(large);
    x.push_back(A::one());
    x.push_back(-large);
    return x;
}

template <class REAL>
void test_verified_regular_cases(const char* tier, const int* exponents, std::size_t exponent_count) {
    using A = vmplapack::Rarith<REAL>;

    REAL eps = exact_epsilon_above_one<REAL>();
    REAL y2 = -A::one() + eps;
    require_vrdot_contains_oracle(vmplapack::gendot::family_a_alternating<REAL>(8, eps), tier);
    require_vrdot_contains_oracle(vmplapack::gendot::family_b_two_term<REAL>(y2), tier);
    require_vrsum_contains_oracle(alternating_sum_case<REAL>(8, eps), tier);
    require_vrsum_contains_oracle(std::vector<REAL>{A::one(), y2}, tier);

    for (std::size_t i = 0; i < exponent_count; ++i) {
        vmplapack::gendot::Rdot_case<REAL> c =
            vmplapack::gendot::family_c_exponent_cancellation<REAL>(exponents[i], A::one());
        require_vrdot_contains_oracle(c, tier);
        require_vrsum_contains_oracle(exponent_sum_case<REAL>(exponents[i]), tier);
    }
}

template <class REAL>
void test_verified_strides(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches verified routines that validate strides but then index as if incx/incy were one.
    std::vector<REAL> xs(static_cast<std::size_t>(13), A::zero());
    std::vector<REAL> ys(static_cast<std::size_t>(10), A::zero());
    xs[0] = A::one();
    xs[3] = -A::one();
    xs[6] = exact_epsilon_above_one<REAL>();
    xs[9] = A::half();
    ys[0] = A::one();
    ys[2] = A::one();
    ys[4] = A::one();
    ys[6] = -A::one();

    vmplapack::Rmidrad<REAL> sum_box = vmplapack::vRsum(static_cast<std::ptrdiff_t>(4), xs.data(), 3);
    std::vector<REAL> sum_y = ones<REAL>(4);
    vmplapack::oracle::Rdot_interval sum_ref =
        vmplapack::oracle::Rdot_oracle(static_cast<std::ptrdiff_t>(4), xs.data(), 3, sum_y.data(), 1);
    require<REAL>(sum_box.status == vmplapack::Rstatus::ok, tier, "strided vRsum returned non-ok");
    require<REAL>(midrad_covers_oracle(sum_box, sum_ref), tier, "strided vRsum did not cover oracle interval");

    vmplapack::Rmidrad<REAL> dot_box = vmplapack::vRdot(static_cast<std::ptrdiff_t>(4), xs.data(), 3, ys.data(), 2);
    vmplapack::oracle::Rdot_interval dot_ref =
        vmplapack::oracle::Rdot_oracle(static_cast<std::ptrdiff_t>(4), xs.data(), 3, ys.data(), 2);
    require<REAL>(dot_box.status == vmplapack::Rstatus::ok, tier, "strided vRdot returned non-ok");
    require<REAL>(midrad_covers_oracle(dot_box, dot_ref), tier, "strided vRdot did not cover oracle interval");
}

template <class REAL>
void test_verified_boundaries_and_statuses(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    // Catches verified empty/invalid input behavior that accidentally follows the accurate UB contract.
    vmplapack::Rmidrad<REAL> empty_sum = vmplapack::vRsum<REAL>(0, nullptr, 1);
    require<REAL>(empty_sum.status == vmplapack::Rstatus::ok && empty_sum.mid == A::zero() && empty_sum.rad == A::zero(),
                  tier,
                  "vRsum n==0 did not return the exact zero enclosure");
    vmplapack::Rmidrad<REAL> empty_dot = vmplapack::vRdot<REAL>(0, nullptr, 1, nullptr, 1);
    require<REAL>(empty_dot.status == vmplapack::Rstatus::ok && empty_dot.mid == A::zero() && empty_dot.rad == A::zero(),
                  tier,
                  "vRdot n==0 did not return the exact zero enclosure");
    require<REAL>(vmplapack::vRsum<REAL>(-1, nullptr, 1).status == vmplapack::Rstatus::invalid_input,
                  tier,
                  "vRsum n<0 was not invalid_input");
    require<REAL>(vmplapack::vRsum<REAL>(1, nullptr, 1).status == vmplapack::Rstatus::invalid_input,
                  tier,
                  "vRsum null pointer was not invalid_input");
    REAL one = A::one();
    require<REAL>(vmplapack::vRsum<REAL>(1, &one, 0).status == vmplapack::Rstatus::invalid_input,
                  tier,
                  "vRsum non-positive stride was not invalid_input");
    require<REAL>(vmplapack::vRdot<REAL>(-1, nullptr, 1, nullptr, 1).status == vmplapack::Rstatus::invalid_input,
                  tier,
                  "vRdot n<0 was not invalid_input");
    require<REAL>(vmplapack::vRdot<REAL>(1, &one, 1, nullptr, 1).status == vmplapack::Rstatus::invalid_input,
                  tier,
                  "vRdot null y pointer was not invalid_input");
    require<REAL>(vmplapack::vRdot<REAL>(1, &one, 0, &one, 1).status == vmplapack::Rstatus::invalid_input,
                  tier,
                  "vRdot non-positive stride was not invalid_input");

    // Catches false finite certificates for NaN/Inf inputs.
    REAL nan_value = quiet_nan_value<REAL>();
    std::vector<REAL> sum_nonfinite{A::one(), nan_value};
    require<REAL>(vmplapack::vRsum(static_cast<std::ptrdiff_t>(sum_nonfinite.size()), sum_nonfinite.data(), 1).status ==
                      vmplapack::Rstatus::non_finite,
                  tier,
                  "vRsum non-finite input was not non_finite");
    REAL inf = A::infinity();
    require<REAL>(vmplapack::vRdot<REAL>(1, &inf, 1, &one, 1).status == vmplapack::Rstatus::non_finite,
                  tier,
                  "vRdot non-finite input was not non_finite");

    // Catches a finite-input overflow reported as ok instead of unbounded.
    REAL large = overflow_large_value<REAL>();
    std::vector<REAL> overflowing_sum{large, large};
    vmplapack::Rmidrad<REAL> sum_overflow =
        vmplapack::vRsum(static_cast<std::ptrdiff_t>(overflowing_sum.size()), overflowing_sum.data(), 1);
    require<REAL>(sum_overflow.status == vmplapack::Rstatus::unbounded, tier, "vRsum overflow was not unbounded");
    REAL two = A::one() + A::one();
    vmplapack::Rmidrad<REAL> dot_overflow = vmplapack::vRdot<REAL>(1, &large, 1, &two, 1);
    require<REAL>(dot_overflow.status == vmplapack::Rstatus::unbounded, tier, "vRdot overflow was not unbounded");
}

template <class REAL>
void require_residual_row_contains_oracle(const REAL* A_data,
                                          std::ptrdiff_t lda,
                                          const REAL* x,
                                          const REAL* b,
                                          std::ptrdiff_t row,
                                          std::ptrdiff_t n,
                                          const vmplapack::Rmidrad<REAL>& box,
                                          const char* tier) {
    // Catches the sign error A*x-b instead of b-A*x by comparing each row with a dot oracle.
    std::vector<REAL> oracle_x;
    std::vector<REAL> oracle_y;
    oracle_x.push_back(b[row]);
    oracle_y.push_back(vmplapack::Rarith<REAL>::one());
    for (std::ptrdiff_t col = 0; col < n; ++col) {
        REAL neg_a = -A_data[row * lda + col];
        oracle_x.push_back(neg_a);
        oracle_y.push_back(x[col]);
    }
    vmplapack::oracle::Rdot_interval ref =
        vmplapack::oracle::Rdot_oracle(static_cast<std::ptrdiff_t>(oracle_x.size()),
                                       oracle_x.data(),
                                       1,
                                       oracle_y.data(),
                                       1);
    require<REAL>(box.status == vmplapack::Rstatus::ok, tier, "vRresidual row returned non-ok");
    require<REAL>(midrad_covers_oracle(box, ref), tier, "vRresidual row did not cover oracle interval");
}

template <class REAL>
void test_residual_regular_and_boundaries(const char* tier) {
    using A = vmplapack::Rarith<REAL>;

    std::vector<REAL> matrix;
    matrix.push_back(A::one() + A::one());
    matrix.push_back(A::one());
    matrix.push_back(-A::one());
    REAL four = A::one() + A::one();
    four = four + four;
    matrix.push_back(four);
    std::vector<REAL> x;
    REAL three = A::one() + A::one();
    three = three + A::one();
    x.push_back(three);
    x.push_back(four);
    std::vector<REAL> b;
    REAL twenty = four;
    twenty = twenty * (four + A::one());
    b.push_back(twenty);
    b.push_back(-A::one());
    std::vector<vmplapack::Rmidrad<REAL>> out(static_cast<std::size_t>(2));

    vmplapack::Rstatus status = vmplapack::vRresidual<REAL>(2, 2, matrix.data(), 2, x.data(), b.data(), out.data());
    require<REAL>(status == vmplapack::Rstatus::ok, tier, "vRresidual regular case returned non-ok");
    require_residual_row_contains_oracle(matrix.data(), 2, x.data(), b.data(), 0, 2, out[0], tier);
    require_residual_row_contains_oracle(matrix.data(), 2, x.data(), b.data(), 1, 2, out[1], tier);

    // Catches m==0 and n==0 boundary handling, including the rule that A and x are not read when n==0.
    require<REAL>(vmplapack::vRresidual<REAL>(0, 2, nullptr, 0, nullptr, nullptr, nullptr) == vmplapack::Rstatus::ok,
                  tier,
                  "vRresidual m==0 did not return ok");
    std::vector<REAL> b_only;
    b_only.push_back(A::one());
    b_only.push_back(-A::one());
    std::vector<vmplapack::Rmidrad<REAL>> n0_out(static_cast<std::size_t>(2));
    require<REAL>(vmplapack::vRresidual<REAL>(2, 0, nullptr, 0, nullptr, b_only.data(), n0_out.data()) ==
                      vmplapack::Rstatus::ok,
                  tier,
                  "vRresidual n==0 did not return ok for finite b");
    require<REAL>(n0_out[0].status == vmplapack::Rstatus::ok && n0_out[0].mid == b_only[0] && n0_out[0].rad == A::zero(),
                  tier,
                  "vRresidual n==0 row 0 is wrong");
    require<REAL>(n0_out[1].status == vmplapack::Rstatus::ok && n0_out[1].mid == b_only[1] && n0_out[1].rad == A::zero(),
                  tier,
                  "vRresidual n==0 row 1 is wrong");

    // Catches invalid residual arguments and verifies that no rows are required for m<0/n<0.
    require<REAL>(vmplapack::vRresidual<REAL>(-1, 0, nullptr, 0, nullptr, nullptr, nullptr) ==
                      vmplapack::Rstatus::invalid_input,
                  tier,
                  "vRresidual m<0 was not invalid_input");
    require<REAL>(vmplapack::vRresidual<REAL>(1, -1, nullptr, 0, nullptr, nullptr, nullptr) ==
                      vmplapack::Rstatus::invalid_input,
                  tier,
                  "vRresidual n<0 was not invalid_input");
    require<REAL>(vmplapack::vRresidual<REAL>(1, 1, nullptr, 1, x.data(), b.data(), out.data()) ==
                      vmplapack::Rstatus::invalid_input,
                  tier,
                  "vRresidual null A was not invalid_input");
    require<REAL>(vmplapack::vRresidual<REAL>(1, 2, matrix.data(), 1, x.data(), b.data(), out.data()) ==
                      vmplapack::Rstatus::invalid_input,
                  tier,
                  "vRresidual lda<n was not invalid_input");

    // Catches residual status aggregation across rows.
    std::vector<REAL> bad_b = b_only;
    bad_b[1] = quiet_nan_value<REAL>();
    require<REAL>(vmplapack::vRresidual<REAL>(2, 0, nullptr, 0, nullptr, bad_b.data(), n0_out.data()) ==
                      vmplapack::Rstatus::non_finite,
                  tier,
                  "vRresidual did not return worst non_finite status");

    // Catches residual finite-input overflow reported as a finite certificate.
    REAL large = overflow_large_value<REAL>();
    REAL two = A::one() + A::one();
    std::vector<REAL> overflow_matrix;
    overflow_matrix.push_back(-large);
    std::vector<REAL> overflow_x;
    overflow_x.push_back(two);
    std::vector<REAL> overflow_b;
    overflow_b.push_back(A::zero());
    std::vector<vmplapack::Rmidrad<REAL>> overflow_out(static_cast<std::size_t>(1));
    require<REAL>(vmplapack::vRresidual<REAL>(1,
                                             1,
                                             overflow_matrix.data(),
                                             1,
                                             overflow_x.data(),
                                             overflow_b.data(),
                                             overflow_out.data()) == vmplapack::Rstatus::unbounded,
                  tier,
                  "vRresidual overflow was not unbounded");
}

template <class REAL>
void test_tier(const char* tier, const int* exponents, std::size_t exponent_count) {
    test_verified_regular_cases<REAL>(tier, exponents, exponent_count);
    test_verified_strides<REAL>(tier);
    test_verified_boundaries_and_statuses<REAL>(tier);
    test_residual_regular_and_boundaries<REAL>(tier);
}

} // namespace

int main() {
    long expected_mpfr_precision = parse_expected_mpfr_precision();
    if (expected_mpfr_precision != 0) {
        mpfrxx::set_default_precision_bits(static_cast<mpfr_prec_t>(expected_mpfr_precision));
    }

    const int float_exponents[] = {24, 30};
    const int double_exponents[] = {53, 60};
    test_tier<float>("float", float_exponents, 2);
    test_tier<double>("double", double_exponents, 2);

    if (expected_mpfr_precision == 53) {
        test_tier<mpfrxx::mpfr_class>("mpfr", double_exponents, 2);
    } else if (expected_mpfr_precision == 512) {
        const int mpfr_w512_exponents[] = {520, 600, 700};
        test_tier<mpfrxx::mpfr_class>("mpfr", mpfr_w512_exponents, 3);
    } else if (expected_mpfr_precision != 0) {
        std::cerr << "unsupported MPFRXX_DEFAULT_PRECISION_BITS=" << expected_mpfr_precision << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
