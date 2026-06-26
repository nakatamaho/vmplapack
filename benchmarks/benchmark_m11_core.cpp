// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include "../tests/Rdot_oracle.h"
#include "../tests/Rgendot.h"

#include <vmplapack/vmplapack.h>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>

#ifndef VMPLAPACK_BENCHMARK_COMPILER_ID
#define VMPLAPACK_BENCHMARK_COMPILER_ID "unknown"
#endif

#ifndef VMPLAPACK_BENCHMARK_COMPILER_VERSION
#define VMPLAPACK_BENCHMARK_COMPILER_VERSION "unknown"
#endif

#ifndef VMPLAPACK_BENCHMARK_BUILD_TYPE
#define VMPLAPACK_BENCHMARK_BUILD_TYPE "unknown"
#endif

#ifndef VMPLAPACK_BENCHMARK_FP_CONTRACT_FLAGS
#define VMPLAPACK_BENCHMARK_FP_CONTRACT_FLAGS "strict-fp"
#endif

#ifndef VMPLAPACK_BENCHMARK_CPU
#define VMPLAPACK_BENCHMARK_CPU "unknown"
#endif

#ifndef VMPLAPACK_BENCHMARK_OS
#define VMPLAPACK_BENCHMARK_OS "unknown"
#endif

#ifndef VMPLAPACK_BENCHMARK_GIT_SHA
#define VMPLAPACK_BENCHMARK_GIT_SHA "unknown"
#endif

namespace {

struct condition_info {
    mpfrxx::mpfr_class cond_oro;
    mpfrxx::mpfr_class cond_sum;
};

struct case_metadata {
    const char* tier;
    const char* generator;
    std::string target_cond;
    std::string seed;
    std::ptrdiff_t n;
    mpfrxx::mpfr_class realized_cond_oro;
    mpfrxx::mpfr_class realized_cond_sum;
    int repetitions;
};

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

mpfrxx::mpfr_class abs_diff(const mpfrxx::mpfr_class& a, const mpfrxx::mpfr_class& b, mpfr_prec_t precision) {
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

std::string format_double(double value) {
    std::ostringstream stream;
    stream << std::setprecision(12) << value;
    return stream.str();
}

std::string format_mpfr(const mpfrxx::mpfr_class& value) {
    if (mpfr_nan_p(value.mpfr_data()) != 0) {
        return "nan";
    }
    if (mpfr_inf_p(value.mpfr_data()) != 0) {
        return (mpfr_sgn(value.mpfr_data()) < 0) ? "-inf" : "+inf";
    }
    std::ostringstream stream;
    stream << std::setprecision(90) << value;
    return stream.str();
}

template <class REAL>
std::string format_real(const REAL& value, mpfr_prec_t precision) {
    mpfrxx::mpfr_class widened = vmplapack::oracle::widen_value(value, precision);
    return format_mpfr(widened);
}

template <class REAL>
mpfrxx::mpfr_class value_error_abs(const REAL& value,
                                   const mpfrxx::mpfr_class& oracle_mid,
                                   mpfr_prec_t precision) {
    mpfrxx::mpfr_class value_mp = vmplapack::oracle::widen_value(value, precision);
    return abs_diff(value_mp, oracle_mid, precision);
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

condition_info compute_conditions(std::ptrdiff_t n,
                                  const mpfrxx::mpfr_class& oracle_mid,
                                  const mpfrxx::mpfr_class& s_up,
                                  mpfr_prec_t precision) {
    condition_info info;
    info.cond_oro = vmplapack::oracle::zero_at(precision);
    info.cond_sum = vmplapack::oracle::zero_at(precision);

    mpfrxx::mpfr_class abs_mid = vmplapack::oracle::abs_at(oracle_mid, precision);
    if (abs_mid == vmplapack::oracle::zero_at(precision)) {
        if (s_up == vmplapack::oracle::zero_at(precision)) {
            return info;
        }
        info.cond_sum = vmplapack::gendot::infinity_at(precision);
        info.cond_oro = vmplapack::gendot::infinity_at(precision);
        return info;
    }

    (void)n;
    mpfrxx::rounding_mode_scope scope(MPFR_RNDU);
    info.cond_sum = s_up / abs_mid;
    mpfrxx::mpfr_class two = mpfrxx::mpfr_class::with_precision(precision, 2.0);
    info.cond_oro = two * info.cond_sum;
    return info;
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

std::int64_t elapsed_ns(std::chrono::steady_clock::time_point start,
                        std::chrono::steady_clock::time_point finish) {
    std::chrono::nanoseconds elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start);
    return static_cast<std::int64_t>(elapsed.count());
}

template <class REAL>
std::int64_t time_naive(const vmplapack::gendot::Rdot_case<REAL>& c, int repetitions, REAL& result) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    for (int r = 0; r < repetitions; ++r) {
        result = naive_dot(c);
    }
    std::chrono::steady_clock::time_point finish = std::chrono::steady_clock::now();
    return elapsed_ns(start, finish);
}

template <class REAL>
std::int64_t time_rdot(const vmplapack::gendot::Rdot_case<REAL>& c, int repetitions, REAL& result) {
    std::ptrdiff_t n = static_cast<std::ptrdiff_t>(c.x.size());
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    for (int r = 0; r < repetitions; ++r) {
        result = vmplapack::Rdot(n, c.x.data(), 1, c.y.data(), 1);
    }
    std::chrono::steady_clock::time_point finish = std::chrono::steady_clock::now();
    return elapsed_ns(start, finish);
}

template <class REAL>
std::int64_t time_vrdot(const vmplapack::gendot::Rdot_case<REAL>& c,
                        int repetitions,
                        vmplapack::Rmidrad<REAL>& result) {
    using A = vmplapack::Rarith<REAL>;

    std::ptrdiff_t n = static_cast<std::ptrdiff_t>(c.x.size());
    result = {A::zero(), A::zero(), vmplapack::Rstatus::ok};
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    for (int r = 0; r < repetitions; ++r) {
        result = vmplapack::vRdot(n, c.x.data(), 1, c.y.data(), 1);
    }
    std::chrono::steady_clock::time_point finish = std::chrono::steady_clock::now();
    return elapsed_ns(start, finish);
}

template <class REAL>
std::int64_t time_vrdot_apriori(const vmplapack::gendot::Rdot_case<REAL>& c,
                                int repetitions,
                                vmplapack::Rmidrad<REAL>& result) {
    using A = vmplapack::Rarith<REAL>;

    std::ptrdiff_t n = static_cast<std::ptrdiff_t>(c.x.size());
    result = {A::zero(), A::zero(), vmplapack::Rstatus::ok};
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    for (int r = 0; r < repetitions; ++r) {
        result = vmplapack::vRdot_apriori(n, c.x.data(), 1, c.y.data(), 1);
    }
    std::chrono::steady_clock::time_point finish = std::chrono::steady_clock::now();
    return elapsed_ns(start, finish);
}

double time_per_item(std::int64_t total_ns, std::ptrdiff_t work_items, int repetitions) {
    double denominator = static_cast<double>(work_items) * static_cast<double>(repetitions);
    if (denominator <= 0.0) {
        return 0.0;
    }
    return static_cast<double>(total_ns) / denominator;
}

void emit_header() {
    std::cout << "schema_version,routine,tier,precision_bits,n,generator,seed,target_cond,"
              << "realized_cond_oro,realized_cond_sum,status,enclosed,elapsed_ns_total,"
              << "time_ns_per_item,work_items,repetitions,statistic,mid_error_abs,mid_error_rel,"
              << "radius,relative_radius,compiler_id,compiler_version,build_type,fp_contract_flags,"
              << "cpu,os,git_sha,rounding_backend" << '\n';
}

void emit_row_prefix(const case_metadata& meta,
                     const char* routine,
                     long precision_bits,
                     const char* status,
                     const char* enclosed,
                     std::int64_t total_ns) {
    std::cout << "m11.1,"
              << routine << ','
              << meta.tier << ','
              << precision_bits << ','
              << meta.n << ','
              << meta.generator << ','
              << meta.seed << ','
              << meta.target_cond << ','
              << format_mpfr(meta.realized_cond_oro) << ','
              << format_mpfr(meta.realized_cond_sum) << ','
              << status << ','
              << enclosed << ','
              << total_ns << ','
              << format_double(time_per_item(total_ns, meta.n, meta.repetitions)) << ','
              << meta.n << ','
              << meta.repetitions << ",mean,";
}

void emit_row_suffix(const char* rounding_backend) {
    std::cout << ','
              << VMPLAPACK_BENCHMARK_COMPILER_ID << ','
              << VMPLAPACK_BENCHMARK_COMPILER_VERSION << ','
              << VMPLAPACK_BENCHMARK_BUILD_TYPE << ','
              << VMPLAPACK_BENCHMARK_FP_CONTRACT_FLAGS << ','
              << VMPLAPACK_BENCHMARK_CPU << ','
              << VMPLAPACK_BENCHMARK_OS << ','
              << VMPLAPACK_BENCHMARK_GIT_SHA << ','
              << rounding_backend << '\n';
}

template <class REAL>
void emit_scalar_row(const case_metadata& meta,
                     const char* routine,
                     const REAL& value,
                     std::int64_t total_ns,
                     const mpfrxx::mpfr_class& oracle_mid,
                     mpfr_prec_t precision) {
    mpfrxx::mpfr_class err_abs = value_error_abs(value, oracle_mid, precision);
    mpfrxx::mpfr_class err_rel = relative_from_abs(err_abs, oracle_mid, precision);
    emit_row_prefix(meta, routine, vmplapack::Rarith<REAL>::precision_bits(), "ok", "null", total_ns);
    std::cout << format_mpfr(err_abs) << ','
              << format_mpfr(err_rel) << ",null,null";
    emit_row_suffix(std::is_same<REAL, mpfrxx::mpfr_class>::value ? "mpfr" : "fenv");
}

template <class REAL>
void emit_verified_row(const case_metadata& meta,
                       const char* routine,
                       const vmplapack::Rmidrad<REAL>& box,
                       std::int64_t total_ns,
                       const vmplapack::oracle::Rdot_interval& ref,
                       const mpfrxx::mpfr_class& oracle_mid,
                       mpfr_prec_t precision) {
    mpfrxx::mpfr_class err_abs = value_error_abs(box.mid, oracle_mid, precision);
    mpfrxx::mpfr_class err_rel = relative_from_abs(err_abs, oracle_mid, precision);
    const char* enclosed = midrad_covers_oracle(box, ref) ? "true" : "false";

    std::string radius = "+inf";
    std::string relative_radius = "+inf";
    if (box.status == vmplapack::Rstatus::ok) {
        mpfrxx::mpfr_class rad_mp = vmplapack::oracle::widen_value(box.rad, precision);
        mpfrxx::mpfr_class rel_rad = relative_from_abs(rad_mp, oracle_mid, precision);
        radius = format_mpfr(rad_mp);
        relative_radius = format_mpfr(rel_rad);
    }

    emit_row_prefix(meta,
                    routine,
                    vmplapack::Rarith<REAL>::precision_bits(),
                    status_name(box.status),
                    enclosed,
                    total_ns);
    std::cout << format_mpfr(err_abs) << ','
              << format_mpfr(err_rel) << ','
              << radius << ','
              << relative_radius;
    emit_row_suffix(std::is_same<REAL, mpfrxx::mpfr_class>::value ? "mpfr" : "fenv");
}

template <class REAL>
void emit_case(const char* tier,
               const char* generator,
               const std::string& target_cond,
               const std::string& seed,
               const vmplapack::gendot::Rdot_case<REAL>& c,
               int repetitions) {
    using A = vmplapack::Rarith<REAL>;

    std::ptrdiff_t n = static_cast<std::ptrdiff_t>(c.x.size());
    vmplapack::oracle::Rdot_interval ref = vmplapack::oracle::Rdot_oracle(n, c.x.data(), 1, c.y.data(), 1);
    mpfr_prec_t precision = ref.precision + static_cast<mpfr_prec_t>(A::precision_bits() + 192L);
    mpfrxx::mpfr_class oracle_mid = interval_midpoint(ref, precision);
    mpfrxx::mpfr_class s_up = vmplapack::oracle::upward_abs_term_sum(n, c.x.data(), 1, c.y.data(), 1, precision);
    condition_info cond = compute_conditions(n, oracle_mid, s_up, precision);

    case_metadata meta;
    meta.tier = tier;
    meta.generator = generator;
    meta.target_cond = target_cond;
    meta.seed = seed;
    meta.n = n;
    meta.realized_cond_oro = cond.cond_oro;
    meta.realized_cond_sum = cond.cond_sum;
    meta.repetitions = repetitions;

    REAL naive = A::zero();
    REAL accurate = A::zero();
    vmplapack::Rmidrad<REAL> directed;
    vmplapack::Rmidrad<REAL> apriori;

    std::int64_t naive_ns = time_naive(c, repetitions, naive);
    std::int64_t rdot_ns = time_rdot(c, repetitions, accurate);
    std::int64_t directed_ns = time_vrdot(c, repetitions, directed);
    std::int64_t apriori_ns = time_vrdot_apriori(c, repetitions, apriori);

    emit_scalar_row(meta, "naive_dot", naive, naive_ns, oracle_mid, precision);
    emit_scalar_row(meta, "Rdot", accurate, rdot_ns, oracle_mid, precision);
    emit_verified_row(meta, "vRdot", directed, directed_ns, ref, oracle_mid, precision);
    emit_verified_row(meta, "vRdot_apriori", apriori, apriori_ns, ref, oracle_mid, precision);
}

template <class REAL>
void emit_random_case(const char* tier,
                      int target_power2,
                      std::uint64_t seed,
                      vmplapack::gendot::Rpermutation permutation,
                      int repetitions) {
    vmplapack::gendot::Rgenerated_dot_case<REAL> generated =
        vmplapack::gendot::randomized_high_condition_power2<REAL>(target_power2, seed, permutation);
    if (generated.status != vmplapack::gendot::Rgenerator_status::ok) {
        std::cerr << tier << ": random benchmark generator did not produce an ok case\n";
        std::exit(EXIT_FAILURE);
    }
    std::string target = format_mpfr(generated.target_cond_oro);
    emit_case(tier, "m8_random_power2", target, std::to_string(static_cast<unsigned long long>(seed)), generated.data, repetitions);
}

template <class REAL>
void emit_tier(const char* tier,
               int family_exponent,
               int smoke_random_target_power2,
               int full_random_target_power2_a,
               int full_random_target_power2_b,
               int adversarial_large_exponent,
               int adversarial_small_exponent,
               int repetitions,
               bool full_mode) {
    using A = vmplapack::Rarith<REAL>;

    vmplapack::gendot::Rdot_case<REAL> family_c =
        vmplapack::gendot::family_c_exponent_cancellation<REAL>(family_exponent, A::one());
    emit_case(tier, "m3_family_c", "null", "null", family_c, repetitions);

    emit_random_case<REAL>(tier,
                           smoke_random_target_power2,
                           11U,
                           vmplapack::gendot::Rpermutation::sorted,
                           repetitions);

    vmplapack::gendot::Rgenerated_dot_case<REAL> huge_small =
        vmplapack::gendot::adversarial_huge_small_mixing<REAL>(adversarial_large_exponent,
                                                               adversarial_small_exponent,
                                                               107U,
                                                               vmplapack::gendot::Rpermutation::shuffled);
    if (huge_small.status != vmplapack::gendot::Rgenerator_status::ok) {
        std::cerr << tier << ": adversarial benchmark generator did not produce an ok case\n";
        std::exit(EXIT_FAILURE);
    }
    emit_case(tier, "m8_adversarial_huge_small", "null", "107", huge_small.data, repetitions);

    if (full_mode) {
        emit_random_case<REAL>(tier,
                               full_random_target_power2_a,
                               29U,
                               vmplapack::gendot::Rpermutation::generated,
                               repetitions);
        emit_random_case<REAL>(tier,
                               full_random_target_power2_b,
                               97U,
                               vmplapack::gendot::Rpermutation::reversed,
                               repetitions);
    }
}

bool parse_full_mode(int argc, char** argv) {
    bool full_mode = false;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--full") == 0) {
            full_mode = true;
        } else if (std::strcmp(argv[i], "--smoke") == 0) {
            full_mode = false;
        } else if (std::strcmp(argv[i], "--help") == 0) {
            std::cout << "usage: benchmark_m11_core [--smoke|--full]\n";
            std::exit(EXIT_SUCCESS);
        } else {
            std::cerr << "unknown benchmark option: " << argv[i] << '\n';
            std::exit(EXIT_FAILURE);
        }
    }
    return full_mode;
}

} // namespace

int main(int argc, char** argv) {
    bool full_mode = parse_full_mode(argc, argv);
    int repetitions = full_mode ? 256 : 32;

    mpfrxx::set_default_precision_bits(512);

    emit_header();
    emit_tier<float>("float", 30, 20, 30, 40, 30, -10, repetitions, full_mode);
    emit_tier<double>("double", 60, 40, 50, 60, 60, -20, repetitions, full_mode);
    emit_tier<mpfrxx::mpfr_class>("mpfr@512", 700, 128, 384, 700, 700, -80, repetitions, full_mode);
    return EXIT_SUCCESS;
}
