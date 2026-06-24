// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include "Rdot_oracle.h"
#include "Rgendot.h"

#include <vmplapack/vmplapack.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

namespace {

/*
 * This test exercises the M8 generator layer in four steps:
 * 1. Build seeded high-condition dot cases at several target ORO condition numbers.
 * 2. Check the generator metadata against the MPFR oracle condition estimate.
 * 3. Run vRdot and vRdot_apriori on the generated cases and require oracle interval inclusion.
 * 4. Check that naive dot and Dot2/Rdot errors fit their expected gamma_n-scaled regimes.
 *
 * Set VMPLAPACK_TEST_VERBOSE=1 when running this executable directly to print the targeted random
 * high-condition cases, measured condition numbers, oracle interval, and verified boxes.
 * Optionally set VMPLAPACK_TEST_VERBOSE_LIMIT=N to print only the first N targeted cases.
 * Optionally set VMPLAPACK_TEST_VERBOSE_TIER=float|double|mpfr to print only one tier.
 */

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

bool verbose_enabled() {
    const char* value = std::getenv("VMPLAPACK_TEST_VERBOSE");
    return value != nullptr && value[0] != '\0' && value[0] != '0';
}

int verbose_case_limit() {
    const char* value = std::getenv("VMPLAPACK_TEST_VERBOSE_LIMIT");
    if (value == nullptr || value[0] == '\0') {
        return -1;
    }

    char* end = nullptr;
    long parsed = std::strtol(value, &end, 10);
    if (end == value || *end != '\0' || parsed < 0) {
        std::cerr << "invalid VMPLAPACK_TEST_VERBOSE_LIMIT=" << value << '\n';
        std::exit(EXIT_FAILURE);
    }
    if (parsed > static_cast<long>(std::numeric_limits<int>::max())) {
        return std::numeric_limits<int>::max();
    }
    return static_cast<int>(parsed);
}

const char* verbose_tier_filter() {
    const char* value = std::getenv("VMPLAPACK_TEST_VERBOSE_TIER");
    if (value == nullptr || value[0] == '\0') {
        return nullptr;
    }
    if (std::strcmp(value, "float") != 0 && std::strcmp(value, "double") != 0 && std::strcmp(value, "mpfr") != 0) {
        std::cerr << "invalid VMPLAPACK_TEST_VERBOSE_TIER=" << value << '\n';
        std::exit(EXIT_FAILURE);
    }
    return value;
}

bool verbose_tier_matches(const char* tier) {
    const char* filter = verbose_tier_filter();
    return filter == nullptr || std::strcmp(filter, tier) == 0;
}

bool verbose_case_should_print(const char* tier) {
    if (!verbose_enabled() || !verbose_tier_matches(tier)) {
        return false;
    }

    static int printed_cases = 0;
    int limit = verbose_case_limit();
    if (limit >= 0 && printed_cases >= limit) {
        return false;
    }
    ++printed_cases;
    return true;
}

const char* permutation_name(vmplapack::gendot::Rpermutation permutation) {
    switch (permutation) {
        case vmplapack::gendot::Rpermutation::generated:
            return "generated";
        case vmplapack::gendot::Rpermutation::sorted:
            return "sorted";
        case vmplapack::gendot::Rpermutation::reversed:
            return "reversed";
        case vmplapack::gendot::Rpermutation::shuffled:
            return "shuffled";
    }
    return "generated";
}

const char* generator_status_name(vmplapack::gendot::Rgenerator_status status) {
    switch (status) {
        case vmplapack::gendot::Rgenerator_status::ok:
            return "ok";
        case vmplapack::gendot::Rgenerator_status::invalid_target:
            return "invalid_target";
        case vmplapack::gendot::Rgenerator_status::unachievable:
            return "unachievable";
    }
    return "unachievable";
}

const char* rstatus_name(vmplapack::Rstatus status) {
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
int output_digits() {
    return std::numeric_limits<REAL>::max_digits10;
}

template <>
int output_digits<mpfrxx::mpfr_class>() {
    return 60;
}

template <class REAL>
void print_generated_case(const char* phase,
                          int target_power2,
                          const vmplapack::gendot::Rgenerated_dot_case<REAL>& generated,
                          bool print_case) {
    if (!print_case) {
        return;
    }

    std::cout << std::setprecision(output_digits<REAL>());
    std::cout << '\n' << phase << '\n';
    std::cout << "  tier = " << generated.tier << '\n';
    if (target_power2 >= 0) {
        std::cout << "  target cond_oro = 2^" << target_power2 << '\n';
    }
    std::cout << "  generator status = " << generator_status_name(generated.status) << '\n';
    std::cout << "  seed = " << generated.seed << '\n';
    std::cout << "  scale = " << generated.scale << '\n';
    std::cout << "  permutation = " << permutation_name(generated.permutation) << '\n';
    std::cout << "  n = " << generated.data.x.size() << '\n';
    std::cout << "  exact = " << generated.data.exact << '\n';
    std::cout << "  target_cond_oro = " << generated.target_cond_oro << '\n';
    std::cout << "  measured_cond_oro = " << generated.measured_cond_oro << '\n';
    std::cout << "  measured_cond_sum = " << generated.measured_cond_sum << '\n';
}

template <class REAL>
void print_midrad_box(const char* label, const vmplapack::Rmidrad<REAL>& box, bool print_case) {
    if (!print_case) {
        return;
    }

    std::cout << std::setprecision(output_digits<REAL>());
    std::cout << "  " << label << '\n';
    std::cout << "    mid = " << box.mid << '\n';
    std::cout << "    rad = " << box.rad << '\n';
    std::cout << "    status = " << rstatus_name(box.status) << '\n';
}

void print_oracle_interval(const vmplapack::oracle::Rdot_interval& ref, bool print_case) {
    if (!print_case) {
        return;
    }

    std::cout << std::setprecision(60);
    std::cout << "  oracle interval" << '\n';
    std::cout << "    lo = " << ref.lo << '\n';
    std::cout << "    hi = " << ref.hi << '\n';
    std::cout << "    precision = " << ref.precision << '\n';
    std::cout << "    refinements = " << ref.refinements << '\n';
}

template <class REAL>
void require_mpfr_case_precision(const vmplapack::gendot::Rdot_case<REAL>&, long) {}

template <>
void require_mpfr_case_precision<mpfrxx::mpfr_class>(const vmplapack::gendot::Rdot_case<mpfrxx::mpfr_class>& c,
                                                     long precision) {
    require<mpfrxx::mpfr_class>(vmplapack::gendot::all_inputs_have_precision(c, precision),
                                "mpfr",
                                "generated inputs are not at W");
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

mpfrxx::mpfr_class log10_value(const mpfrxx::mpfr_class& value, mpfr_prec_t precision) {
    mpfrxx::mpfr_class widened = vmplapack::oracle::widen_mpfr(value, precision);
    mpfrxx::mpfr_class result = mpfrxx::mpfr_class::with_precision(precision);
    mpfr_log10(result.mpfr_data(), widened.mpfr_data(), MPFR_RNDN);
    return result;
}

mpfrxx::mpfr_class abs_diff(const mpfrxx::mpfr_class& a, const mpfrxx::mpfr_class& b, mpfr_prec_t precision) {
    mpfrxx::mpfr_class aa = vmplapack::oracle::widen_mpfr(a, precision);
    mpfrxx::mpfr_class bb = vmplapack::oracle::widen_mpfr(b, precision);
    mpfrxx::mpfr_class result = mpfrxx::mpfr_class::with_precision(precision);
    if (aa >= bb) {
        result = aa - bb;
    } else {
        result = bb - aa;
    }
    return result;
}

bool log10_close(const mpfrxx::mpfr_class& measured,
                 const mpfrxx::mpfr_class& target,
                 double tolerance_decades) {
    if (mpfr_number_p(measured.mpfr_data()) == 0 || mpfr_number_p(target.mpfr_data()) == 0) {
        return false;
    }
    mpfr_prec_t precision = std::max(mpfr_get_prec(measured.mpfr_data()), mpfr_get_prec(target.mpfr_data())) + 128;
    mpfrxx::mpfr_class log_measured = log10_value(measured, precision);
    mpfrxx::mpfr_class log_target = log10_value(target, precision);
    mpfrxx::mpfr_class diff = abs_diff(log_measured, log_target, precision);
    mpfrxx::mpfr_class tolerance = mpfrxx::mpfr_class::with_precision(precision, tolerance_decades);
    return diff <= tolerance;
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
mpfrxx::mpfr_class relative_error(const REAL& got,
                                  const vmplapack::oracle::Rdot_interval& ref,
                                  mpfr_prec_t precision) {
    mpfrxx::mpfr_class got_mp = vmplapack::oracle::widen_value(got, precision);
    mpfrxx::mpfr_class ref_mid = interval_midpoint(ref, precision);
    mpfrxx::mpfr_class diff = abs_diff(got_mp, ref_mid, precision);
    mpfrxx::mpfr_class abs_ref = vmplapack::oracle::abs_at(ref_mid, precision);
    if (abs_ref == vmplapack::oracle::zero_at(precision)) {
        return diff;
    }
    mpfrxx::rounding_mode_scope scope(MPFR_RNDU);
    mpfrxx::mpfr_class rel = diff / abs_ref;
    return rel;
}

template <class REAL>
void require_condition_metadata(const vmplapack::gendot::Rgenerated_dot_case<REAL>& generated,
                                const char* tier) {
    // Phase 1: verify generator metadata before any kernel is tested.
    // Catches the factor-of-two mismatch between ORO cond=2*S/|dot| and the sum convention S/|dot|.
    require<REAL>(generated.status == vmplapack::gendot::Rgenerator_status::ok,
                  tier,
                  "generator did not produce an ok case");
    require<REAL>(log10_close(generated.measured_cond_oro, generated.target_cond_oro, 0.25),
                  tier,
                  "measured ORO condition did not match target within default log tolerance");
    mpfr_prec_t precision = std::max(mpfr_get_prec(generated.measured_cond_oro.mpfr_data()),
                                     mpfr_get_prec(generated.measured_cond_sum.mpfr_data())) + 64;
    mpfrxx::mpfr_class two = mpfrxx::mpfr_class::with_precision(precision, 2.0);
    mpfrxx::mpfr_class twice_sum = two * vmplapack::oracle::widen_mpfr(generated.measured_cond_sum, precision);
    mpfrxx::mpfr_class oro = vmplapack::oracle::widen_mpfr(generated.measured_cond_oro, precision);
    require<REAL>(twice_sum == oro, tier, "measured cond_oro was not exactly twice cond_sum");
}

template <class REAL>
void require_verified_inclusion(const vmplapack::gendot::Rdot_case<REAL>& c, const char* tier, bool print_case) {
    // Phase 2: use the generated data as a verified-dot regression case.
    // Catches a generated ordering that breaks directed-rounding verified inclusion.
    vmplapack::Rmidrad<REAL> box = vmplapack::vRdot(static_cast<std::ptrdiff_t>(c.x.size()),
                                                   c.x.data(),
                                                   1,
                                                   c.y.data(),
                                                   1);
    vmplapack::oracle::Rdot_interval ref = vmplapack::oracle::Rdot_oracle(static_cast<std::ptrdiff_t>(c.x.size()),
                                                                         c.x.data(),
                                                                         1,
                                                                         c.y.data(),
                                                                         1);
    print_oracle_interval(ref, print_case);
    print_midrad_box("vRdot", box, print_case);
    require<REAL>(box.status == vmplapack::Rstatus::ok, tier, "vRdot generated case returned non-ok");
    require<REAL>(midrad_covers_oracle(box, ref), tier, "vRdot did not cover generated oracle interval");
}

template <class REAL>
void require_apriori_inclusion(const vmplapack::gendot::Rdot_case<REAL>& c, const char* tier, bool print_case) {
    // Phase 3: reuse the same generated data against the M7 a-priori certificate.
    // Catches M8-generated high-condition cases where the M7 a-priori certificate misses the oracle interval.
    std::ptrdiff_t n = static_cast<std::ptrdiff_t>(c.x.size());
    vmplapack::Rmidrad<REAL> box = vmplapack::vRdot_apriori(n, c.x.data(), 1, c.y.data(), 1);
    vmplapack::oracle::Rdot_interval ref = vmplapack::oracle::Rdot_oracle(n, c.x.data(), 1, c.y.data(), 1);
    print_midrad_box("vRdot_apriori", box, print_case);
    require<REAL>(box.status == vmplapack::Rstatus::ok, tier, "vRdot_apriori generated case returned non-ok");
    require<REAL>(box.mid == vmplapack::Rdot(n, c.x.data(), 1, c.y.data(), 1),
                  tier,
                  "vRdot_apriori generated case midpoint is not Rdot");
    require<REAL>(midrad_covers_oracle(box, ref),
                  tier,
                  "vRdot_apriori did not cover generated oracle interval");
}

template <class REAL>
void require_error_scaling(const vmplapack::gendot::Rgenerated_dot_case<REAL>& generated,
                           const char* tier) {
    // Phase 4: check that the generated condition number drives the expected error scale.
    // Catches tests that use the old cond~1/u threshold instead of the gamma_n-scaled error shapes.
    const vmplapack::gendot::Rdot_case<REAL>& c = generated.data;
    std::ptrdiff_t n = static_cast<std::ptrdiff_t>(c.x.size());
    vmplapack::oracle::Rdot_interval ref = vmplapack::oracle::Rdot_oracle(n, c.x.data(), 1, c.y.data(), 1);
    mpfr_prec_t precision = ref.precision + static_cast<mpfr_prec_t>(vmplapack::Rarith<REAL>::precision_bits() + 192L);

    REAL naive = naive_dot(c);
    REAL accurate = vmplapack::Rdot(n, c.x.data(), 1, c.y.data(), 1);
    mpfrxx::mpfr_class naive_rel = relative_error(naive, ref, precision);
    mpfrxx::mpfr_class accurate_rel = relative_error(accurate, ref, precision);

    mpfrxx::mpfr_class u = vmplapack::oracle::power_of_two(precision,
                                                           -static_cast<mpfr_exp_t>(vmplapack::Rarith<REAL>::precision_bits()));
    mpfrxx::mpfr_class n_mp = mpfrxx::mpfr_class::with_precision(precision);
    mpfr_set_sj(n_mp.mpfr_data(), static_cast<intmax_t>(n), MPFR_RNDN);
    mpfrxx::mpfr_class one = vmplapack::oracle::one_at(precision);
    mpfrxx::mpfr_class gamma = mpfrxx::mpfr_class::with_precision(precision);
    {
        mpfrxx::rounding_mode_scope scope(MPFR_RNDU);
        mpfrxx::mpfr_class nu = n_mp * u;
        mpfrxx::mpfr_class den = one - nu;
        gamma = nu / den;
    }

    mpfrxx::rounding_mode_scope scope(MPFR_RNDU);
    mpfrxx::mpfr_class naive_scale = gamma * vmplapack::oracle::widen_mpfr(generated.measured_cond_sum, precision);
    mpfrxx::mpfr_class gamma_sq = gamma * gamma;
    mpfrxx::mpfr_class half = mpfrxx::mpfr_class::with_precision(precision, 0.5);
    mpfrxx::mpfr_class dot2_term = gamma_sq * vmplapack::oracle::widen_mpfr(generated.measured_cond_oro, precision);
    mpfrxx::mpfr_class dot2_scale = u + (half * dot2_term);
    mpfrxx::mpfr_class naive_slack = mpfrxx::mpfr_class::with_precision(precision, 8.0) * naive_scale;
    mpfrxx::mpfr_class dot2_slack = mpfrxx::mpfr_class::with_precision(precision, 32.0) * dot2_scale;

    require<REAL>(naive_rel <= naive_slack, tier, "naive error exceeded first-order gamma_n*cond_sum scale");
    require<REAL>(accurate_rel <= dot2_slack, tier, "Rdot error exceeded Dot2 gamma_n^2*cond_oro scale");
}

template <class REAL>
void require_zero_dot_infinite_condition(const char* tier, int scale) {
    // Catches exact-zero adversarial cases being assigned a finite condition number.
    vmplapack::gendot::Rgenerated_dot_case<REAL> zero_case =
        vmplapack::gendot::adversarial_exact_cancellation<REAL>(scale,
                                                                7U,
                                                                vmplapack::gendot::Rpermutation::generated);
    require<REAL>(zero_case.status == vmplapack::gendot::Rgenerator_status::ok,
                  tier,
                  "exact-cancellation generator returned non-ok");
    require<REAL>(mpfr_inf_p(zero_case.measured_cond_oro.mpfr_data()) != 0,
                  tier,
                  "exact nontrivial zero dot did not report infinite ORO condition");
}

template <class REAL>
void test_targeted_random_cases(const char* tier,
                                const int* target_powers,
                                std::size_t target_count,
                                long mpfr_precision) {
    const std::uint64_t seeds[] = {11U, 29U, 97U};
    const vmplapack::gendot::Rpermutation permutations[] = {
        vmplapack::gendot::Rpermutation::generated,
        vmplapack::gendot::Rpermutation::sorted,
        vmplapack::gendot::Rpermutation::reversed,
        vmplapack::gendot::Rpermutation::shuffled};

    for (std::size_t target_index = 0; target_index < target_count; ++target_index) {
        int target_power2 = target_powers[target_index];
        for (std::size_t i = 0; i < 3U; ++i) {
            for (std::size_t j = 0; j < 4U; ++j) {
                vmplapack::gendot::Rgenerated_dot_case<REAL> generated =
                    vmplapack::gendot::randomized_high_condition_power2<REAL>(target_power2,
                                                                              seeds[i],
                                                                              permutations[j]);
                bool print_case = verbose_case_should_print(tier);
                print_generated_case("targeted random high-condition dot", target_power2, generated, print_case);
                require_condition_metadata(generated, tier);
                require_mpfr_case_precision(generated.data, mpfr_precision);
                require_verified_inclusion(generated.data, tier, print_case);
                require_apriori_inclusion(generated.data, tier, print_case);
                require_error_scaling(generated, tier);
            }
        }
    }
}

template <class REAL>
void test_adversarial_families(const char* tier, int scale, int small_exponent, long mpfr_precision) {
    const vmplapack::gendot::Rpermutation permutations[] = {
        vmplapack::gendot::Rpermutation::generated,
        vmplapack::gendot::Rpermutation::sorted,
        vmplapack::gendot::Rpermutation::reversed,
        vmplapack::gendot::Rpermutation::shuffled};

    for (std::size_t i = 0; i < 4U; ++i) {
        vmplapack::gendot::Rgenerated_dot_case<REAL> heavy =
            vmplapack::gendot::adversarial_heavy_cancellation<REAL>(scale, 101U, permutations[i]);
        vmplapack::gendot::Rgenerated_dot_case<REAL> alternating =
            vmplapack::gendot::adversarial_alternating_signs<REAL>(8, small_exponent, 103U, permutations[i]);
        vmplapack::gendot::Rgenerated_dot_case<REAL> huge_small =
            vmplapack::gendot::adversarial_huge_small_mixing<REAL>(scale, small_exponent, 107U, permutations[i]);
        require<REAL>(heavy.status == vmplapack::gendot::Rgenerator_status::ok, tier, "heavy generator failed");
        require<REAL>(alternating.status == vmplapack::gendot::Rgenerator_status::ok, tier, "alternating generator failed");
        require<REAL>(huge_small.status == vmplapack::gendot::Rgenerator_status::ok, tier, "huge/small generator failed");
        require_mpfr_case_precision(heavy.data, mpfr_precision);
        require_mpfr_case_precision(alternating.data, mpfr_precision);
        require_mpfr_case_precision(huge_small.data, mpfr_precision);
        require_verified_inclusion(heavy.data, tier, false);
        require_verified_inclusion(alternating.data, tier, false);
        require_verified_inclusion(huge_small.data, tier, false);
    }
}

template <class REAL>
void test_invalid_and_unachievable(const char* tier) {
    // Catches finite target_cond<2 being accepted under the ORO condition-number convention.
    vmplapack::gendot::Rgenerated_dot_case<REAL> invalid =
        vmplapack::gendot::randomized_high_condition_power2<REAL>(0,
                                                                  5U,
                                                                  vmplapack::gendot::Rpermutation::generated);
    require<REAL>(invalid.status == vmplapack::gendot::Rgenerator_status::invalid_target,
                  tier,
                  "target condition below 2 was not invalid");

    // Catches impossible exponent ranges being silently rounded into a misleading finite target.
    int huge_target = vmplapack::gendot::max_power_exponent<REAL>() + 64;
    vmplapack::gendot::Rgenerated_dot_case<REAL> unachievable =
        vmplapack::gendot::randomized_high_condition_power2<REAL>(huge_target,
                                                                  5U,
                                                                  vmplapack::gendot::Rpermutation::generated);
    require<REAL>(unachievable.status == vmplapack::gendot::Rgenerator_status::unachievable,
                  tier,
                  "unachievable target did not report unachievable");
}

template <class REAL>
void test_tier(const char* tier,
               const int* target_powers,
               std::size_t target_count,
               int scale,
               int small_exponent,
               long mpfr_precision) {
    test_targeted_random_cases<REAL>(tier, target_powers, target_count, mpfr_precision);
    test_adversarial_families<REAL>(tier, scale, small_exponent, mpfr_precision);
    require_zero_dot_infinite_condition<REAL>(tier, scale);
    test_invalid_and_unachievable<REAL>(tier);
}

} // namespace

int main() {
    long expected_mpfr_precision = parse_expected_mpfr_precision();
    if (expected_mpfr_precision != 0) {
        mpfrxx::set_default_precision_bits(static_cast<mpfr_prec_t>(expected_mpfr_precision));
    }

    // Each tier gets low, medium, and high finite ORO-condition targets.
    const int float_targets[] = {10, 20, 30};
    const int double_targets[] = {20, 40, 60};
    test_tier<float>("float", float_targets, 3, 30, -10, 0);
    test_tier<double>("double", double_targets, 3, 60, -20, 0);

    if (expected_mpfr_precision == 53) {
        test_tier<mpfrxx::mpfr_class>("mpfr", double_targets, 3, 60, -20, expected_mpfr_precision);
    } else if (expected_mpfr_precision == 512) {
        const int mpfr_w512_targets[] = {128, 384, 700};
        test_tier<mpfrxx::mpfr_class>("mpfr", mpfr_w512_targets, 3, 700, -80, expected_mpfr_precision);
    } else if (expected_mpfr_precision != 0) {
        std::cerr << "unsupported MPFRXX_DEFAULT_PRECISION_BITS=" << expected_mpfr_precision << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
