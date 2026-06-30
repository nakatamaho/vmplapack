// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include <vmplapack/vmplapack_arith.h>
#include <vmplapack/vmplapack_eft.h>
#include <vmplapack/vmplapack_utils.h>

namespace vmplapack {

template <class REAL>
REAL Rsum(std::ptrdiff_t n, const REAL* x, std::ptrdiff_t incx = 1);

template <class REAL>
REAL Rdot(std::ptrdiff_t n,
          const REAL* x,
          std::ptrdiff_t incx,
          const REAL* y,
          std::ptrdiff_t incy);

template <class REAL>
REAL Rsum_abs_dot_upper(std::ptrdiff_t n,
                        const REAL* x,
                        std::ptrdiff_t incx,
                        const REAL* y,
                        std::ptrdiff_t incy);

template <class REAL>
Rmidrad<REAL> vRsum(std::ptrdiff_t n, const REAL* x, std::ptrdiff_t incx = 1);

template <class REAL>
Rmidrad<REAL> vRdot(std::ptrdiff_t n,
                    const REAL* x,
                    std::ptrdiff_t incx,
                    const REAL* y,
                    std::ptrdiff_t incy);

template <class REAL>
Rmidrad<REAL> vRdot_apriori(std::ptrdiff_t n,
                            const REAL* x,
                            std::ptrdiff_t incx,
                            const REAL* y,
                            std::ptrdiff_t incy);

template <class REAL>
Rstatus vRresidual(std::ptrdiff_t m,
                   std::ptrdiff_t n,
                   const REAL* A,
                   std::ptrdiff_t lda,
                   const REAL* x,
                   const REAL* b,
                   Rmidrad<REAL>* out);

template <class REAL>
VerificationStatus vRgemv_point(std::ptrdiff_t m,
                                std::ptrdiff_t n,
                                const REAL* A,
                                std::ptrdiff_t lda,
                                const REAL* x,
                                std::ptrdiff_t incx,
                                Rmidrad<REAL>* out);

template <class REAL>
VerificationStatus vRgemm_point(std::ptrdiff_t m,
                                std::ptrdiff_t n,
                                std::ptrdiff_t k,
                                const REAL* A,
                                std::ptrdiff_t lda,
                                const REAL* B,
                                std::ptrdiff_t ldb,
                                Rmidrad<REAL>* C,
                                std::ptrdiff_t ldc);

template <class REAL>
VerificationStatus vRgesv(std::ptrdiff_t n,
                          const REAL* A,
                          std::ptrdiff_t lda,
                          const REAL* b,
                          std::ptrdiff_t incb,
                          Rmidrad<REAL>* x,
                          std::ptrdiff_t incx);

template <class REAL>
VerificationStatus vRgesv(std::ptrdiff_t n,
                          std::ptrdiff_t nrhs,
                          const REAL* A,
                          std::ptrdiff_t lda,
                          const REAL* B,
                          std::ptrdiff_t ldb,
                          Rmidrad<REAL>* X,
                          std::ptrdiff_t ldx);

template <class REAL>
VerificationStatus vRgeinv(std::ptrdiff_t n,
                           const REAL* A,
                           std::ptrdiff_t lda,
                           Rmidrad<REAL>* Ainv,
                           std::ptrdiff_t ldinv);

template <class REAL>
struct Rcond_bounds {
    REAL normA_upper;
    REAL normAinv_upper;
    REAL kappa_upper;
    REAL rcond_lower;
    Rstatus status;
};

template <class REAL>
Rcond_bounds<REAL> vRgecon(std::ptrdiff_t n, const REAL* A, std::ptrdiff_t lda);

template <class REAL>
REAL Rsum(std::ptrdiff_t n, const REAL* x, std::ptrdiff_t incx) {
    using A = Rarith<REAL>;

    if (n <= 0) {
        return A::zero();
    }
    assert(x != nullptr);
    assert(incx >= 1);

    REAL s = x[0];
    REAL c = A::zero();
    std::ptrdiff_t ix = incx;
    for (std::ptrdiff_t i = 1; i < n; ++i) {
        REAL e = A::zero();
        Rtwosum(s, x[ix], s, e);
        REAL cn = c + e;
        c = cn;
        ix += incx;
    }
    REAL result = s + c;
    return result;
}

template <class REAL>
REAL Rdot(std::ptrdiff_t n,
          const REAL* x,
          std::ptrdiff_t incx,
          const REAL* y,
          std::ptrdiff_t incy) {
    using A = Rarith<REAL>;

    if (n <= 0) {
        return A::zero();
    }
    assert(x != nullptr);
    assert(y != nullptr);
    assert(incx >= 1);
    assert(incy >= 1);

    REAL p = A::zero();
    REAL s = A::zero();
    Rtwoproduct(x[0], y[0], p, s);

    std::ptrdiff_t ix = incx;
    std::ptrdiff_t iy = incy;
    for (std::ptrdiff_t i = 1; i < n; ++i) {
        REAL h = A::zero();
        REAL r = A::zero();
        Rtwoproduct(x[ix], y[iy], h, r);
        REAL q = A::zero();
        Rtwosum(p, h, p, q);
        REAL t = q + r;
        REAL sn = s + t;
        s = sn;
        ix += incx;
        iy += incy;
    }

    REAL result = p + s;
    return result;
}

template <class REAL>
REAL Rsum_abs_dot_upper(std::ptrdiff_t n,
                        const REAL* x,
                        std::ptrdiff_t incx,
                        const REAL* y,
                        std::ptrdiff_t incy) {
    using A = Rarith<REAL>;

    if (n <= 0) {
        return A::zero();
    }
    assert(x != nullptr);
    assert(y != nullptr);
    assert(incx >= 1);
    assert(incy >= 1);

    typename A::round_up scope;
    REAL acc = A::zero();
    for (std::ptrdiff_t i = 0; i < n; ++i) {
        REAL ax = A::abs(x[i * incx]);
        REAL ay = A::abs(y[i * incy]);
        REAL prod = ax * ay;
        REAL next = acc + prod;
        acc = next;
    }
    return acc;
}

namespace detail {

template <class REAL>
Rmidrad<REAL> invalid_midrad() {
    using A = Rarith<REAL>;
    return {A::zero(), A::infinity(), Rstatus::invalid_input};
}

template <class REAL>
Rmidrad<REAL> non_finite_midrad() {
    using A = Rarith<REAL>;
    return {A::zero(), A::infinity(), Rstatus::non_finite};
}

inline Rstatus worst_status(Rstatus a, Rstatus b) {
    return (Rstatus_rank(a) >= Rstatus_rank(b)) ? a : b;
}

inline VerificationStatus verification_status_from_box(Rstatus status) {
    if (status == Rstatus::ok) {
        return VerificationStatus::Verified;
    }
    if (status == Rstatus::invalid_input) {
        return VerificationStatus::InvalidInput;
    }
    return VerificationStatus::Unverified;
}

inline VerificationStatus worst_verification_status(VerificationStatus a, VerificationStatus b) {
    return (VerificationStatus_rank(a) >= VerificationStatus_rank(b)) ? a : b;
}

template <class REAL>
Rmidrad<REAL> vRdot_apriori_with_reference_fallback(std::ptrdiff_t n,
                                                    const REAL* x,
                                                    std::ptrdiff_t incx,
                                                    const REAL* y,
                                                    std::ptrdiff_t incy) {
    Rmidrad<REAL> fast = vRdot_apriori(n, x, incx, y, incy);
    if (fast.status != Rstatus::unbounded) {
        return fast;
    }

    Rmidrad<REAL> reference = vRdot(n, x, incx, y, incy);
    if (Rstatus_rank(reference.status) < Rstatus_rank(fast.status)) {
        return reference;
    }
    return fast;
}

inline bool apriori_length_gate(std::ptrdiff_t n, long precision_bits) {
    if (n < 0 || precision_bits <= 0) {
        return false;
    }
    long length_digits = static_cast<long>(std::numeric_limits<std::ptrdiff_t>::digits);
    if (precision_bits >= length_digits) {
        return true;
    }
    std::ptrdiff_t limit = std::ptrdiff_t{1} << static_cast<int>(precision_bits);
    return n < limit;
}

template <class REAL>
REAL length_as_real_up(std::ptrdiff_t n) {
    return static_cast<REAL>(n);
}

#ifdef VMPLAPACK_ENABLE_MPFR
template <>
inline mpfrxx::mpfr_class length_as_real_up<mpfrxx::mpfr_class>(std::ptrdiff_t n) {
    mpfrxx::mpfr_class value =
        mpfrxx::mpfr_class::with_precision(static_cast<mpfr_prec_t>(Rarith<mpfrxx::mpfr_class>::precision_bits()));
    mpfr_set_sj(value.mpfr_data(), static_cast<intmax_t>(n), MPFR_RNDU);
    return value;
}
#endif

template <class REAL>
void clear_underflow_flag() {}

template <class REAL>
bool underflow_flag_raised() {
    return false;
}

#ifdef VMPLAPACK_ENABLE_MPFR
template <>
inline void clear_underflow_flag<mpfrxx::mpfr_class>() {
    mpfr_clear_underflow();
}

template <>
inline bool underflow_flag_raised<mpfrxx::mpfr_class>() {
    return mpfr_underflow_p() != 0;
}
#endif

template <class REAL>
Rmidrad<REAL> unbounded_midrad_with_mid(REAL mid) {
    using A = Rarith<REAL>;
    REAL m = A::is_finite(mid) ? mid : A::zero();
    return {m, A::infinity(), Rstatus::unbounded};
}

template <class REAL>
REAL five() {
    using A = Rarith<REAL>;
    REAL two = A::one() + A::one();
    REAL four = two + two;
    REAL value = four + A::one();
    return value;
}

template <class REAL>
REAL Rdot_apriori_radius(std::ptrdiff_t n, REAL mid, REAL S_up) {
    using A = Rarith<REAL>;

    // ORO Prop. 5.5 Dot2 a-priori estimate, not Algorithm 5.8 Dot2Err.
    REAL u = A::unit_roundoff();
    REAL one = A::one();

    REAL n_real = A::zero();
    REAL nu = A::zero();
    {
        typename A::round_up scope;
        n_real = length_as_real_up<REAL>(n);
        REAL nu_product = n_real * u;
        nu = nu_product;
    }

    REAL gamma_den = A::zero();
    {
        typename A::round_down scope;
        REAL den = one - nu;
        gamma_den = den;
    }
    if (gamma_den <= A::zero()) {
        return A::infinity();
    }

    REAL numerator = A::zero();
    {
        typename A::round_up scope;
        REAL gamma = nu / gamma_den;
        REAL gamma_sq = gamma * gamma;
        REAL gamma_term = gamma_sq * S_up;
        REAL abs_mid = A::abs(mid);
        REAL mid_term = u * abs_mid;
        REAL underflow_n = n_real * A::eta();
        REAL underflow_term = underflow_n * five<REAL>();
        REAL partial = mid_term + gamma_term;
        numerator = partial + underflow_term;
    }

    REAL final_den = A::zero();
    {
        typename A::round_down scope;
        REAL den = one - u;
        final_den = den;
    }
    if (final_den <= A::zero()) {
        return A::infinity();
    }

    typename A::round_up scope;
    REAL rad = numerator / final_den;
    return rad;
}

template <class LEFT, class RIGHT>
bool storage_ranges_overlap(const LEFT* left,
                            std::ptrdiff_t left_last_index,
                            const RIGHT* right,
                            std::ptrdiff_t right_last_index) {
    if (left == nullptr || right == nullptr || left_last_index < 0 || right_last_index < 0) {
        return false;
    }
    std::uintptr_t left_begin = reinterpret_cast<std::uintptr_t>(left);
    std::uintptr_t right_begin = reinterpret_cast<std::uintptr_t>(right);
    std::uintptr_t left_size = static_cast<std::uintptr_t>(left_last_index + 1) * sizeof(LEFT);
    std::uintptr_t right_size = static_cast<std::uintptr_t>(right_last_index + 1) * sizeof(RIGHT);
    std::uintptr_t left_end = left_begin + left_size;
    std::uintptr_t right_end = right_begin + right_size;
    return left_begin < right_end && right_begin < left_end;
}

inline std::ptrdiff_t matrix_last_index(std::ptrdiff_t rows, std::ptrdiff_t cols, std::ptrdiff_t ld) {
    if (rows <= 0 || cols <= 0) {
        return -1;
    }
    return (rows - 1) * ld + (cols - 1);
}

inline std::ptrdiff_t strided_last_index(std::ptrdiff_t n, std::ptrdiff_t inc) {
    if (n <= 0) {
        return -1;
    }
    return (n - 1) * inc;
}

template <class REAL>
Rmidrad<REAL> status_midrad(Rstatus status) {
    using A = Rarith<REAL>;
    if (status == Rstatus::ok) {
        return {A::zero(), A::zero(), Rstatus::ok};
    }
    return {A::zero(), A::infinity(), status};
}

template <class REAL>
void fill_strided_boxes(std::ptrdiff_t n, Rmidrad<REAL>* out, std::ptrdiff_t inc, Rstatus status) {
    Rmidrad<REAL> box = status_midrad<REAL>(status);
    for (std::ptrdiff_t i = 0; i < n; ++i) {
        out[i * inc] = box;
    }
}

template <class REAL>
void fill_matrix_boxes(std::ptrdiff_t n, std::ptrdiff_t nrhs, Rmidrad<REAL>* out, std::ptrdiff_t ld, Rstatus status) {
    Rmidrad<REAL> box = status_midrad<REAL>(status);
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        for (std::ptrdiff_t rhs = 0; rhs < nrhs; ++rhs) {
            out[row * ld + rhs] = box;
        }
    }
}

template <class REAL>
bool finite_square_matrix(std::ptrdiff_t n, const REAL* A_data, std::ptrdiff_t lda) {
    using A = Rarith<REAL>;
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            if (!A::is_finite(A_data[row * lda + col])) {
                return false;
            }
        }
    }
    return true;
}

template <class REAL>
Rcond_bounds<REAL> invalid_cond_bounds() {
    using A = Rarith<REAL>;
    return {A::zero(), A::zero(), A::zero(), A::zero(), Rstatus::invalid_input};
}

template <class REAL>
Rcond_bounds<REAL> non_finite_cond_bounds() {
    using A = Rarith<REAL>;
    return {A::infinity(), A::infinity(), A::infinity(), A::zero(), Rstatus::non_finite};
}

template <class REAL>
Rcond_bounds<REAL> unbounded_cond_bounds(REAL normA_upper) {
    using A = Rarith<REAL>;
    REAL normA = A::is_finite(normA_upper) ? normA_upper : A::infinity();
    return {normA, A::infinity(), A::infinity(), A::zero(), Rstatus::unbounded};
}

template <class REAL>
Rcond_bounds<REAL> empty_cond_bounds() {
    using A = Rarith<REAL>;
    return {A::zero(), A::zero(), A::one(), A::one(), Rstatus::ok};
}

template <class REAL>
bool infinity_norm_point_upper(std::ptrdiff_t n, const REAL* A_data, std::ptrdiff_t lda, REAL& norm_upper) {
    using A = Rarith<REAL>;

    norm_upper = A::zero();
    typename A::round_up scope;
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        REAL row_sum = A::zero();
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            REAL term = A::abs(A_data[row * lda + col]);
            REAL next = row_sum + term;
            row_sum = next;
        }
        if (!A::is_finite(row_sum)) {
            return false;
        }
        if (row_sum > norm_upper) {
            norm_upper = row_sum;
        }
    }
    return A::is_finite(norm_upper);
}

template <class REAL>
Rstatus worst_midrad_matrix_status(std::ptrdiff_t n, const Rmidrad<REAL>* boxes, std::ptrdiff_t ld) {
    Rstatus worst = Rstatus::ok;
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            worst = worst_status(worst, boxes[row * ld + col].status);
        }
    }
    return worst;
}

template <class REAL>
bool infinity_norm_midrad_upper(std::ptrdiff_t n, const Rmidrad<REAL>* boxes, std::ptrdiff_t ld, REAL& norm_upper) {
    using A = Rarith<REAL>;

    norm_upper = A::zero();
    typename A::round_up scope;
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        REAL row_sum = A::zero();
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            const Rmidrad<REAL>& box = boxes[row * ld + col];
            if (box.status != Rstatus::ok || box.rad < A::zero() || !A::is_finite(box.mid) || !A::is_finite(box.rad)) {
                return false;
            }
            REAL abs_mid = A::abs(box.mid);
            REAL term = abs_mid + box.rad;
            REAL next = row_sum + term;
            row_sum = next;
        }
        if (!A::is_finite(row_sum)) {
            return false;
        }
        if (row_sum > norm_upper) {
            norm_upper = row_sum;
        }
    }
    return A::is_finite(norm_upper);
}

template <class REAL>
bool finite_rhs_matrix(std::ptrdiff_t n, std::ptrdiff_t nrhs, const REAL* B_data, std::ptrdiff_t ldb) {
    using A = Rarith<REAL>;
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        for (std::ptrdiff_t rhs = 0; rhs < nrhs; ++rhs) {
            if (!A::is_finite(B_data[row * ldb + rhs])) {
                return false;
            }
        }
    }
    return true;
}

template <class REAL>
bool finite_contiguous_matrix(std::ptrdiff_t rows, std::ptrdiff_t cols, const std::vector<REAL>& data) {
    using A = Rarith<REAL>;
    for (std::ptrdiff_t row = 0; row < rows; ++row) {
        for (std::ptrdiff_t col = 0; col < cols; ++col) {
            if (!A::is_finite(data[static_cast<std::size_t>(row * cols + col)])) {
                return false;
            }
        }
    }
    return true;
}

template <class REAL>
bool compute_inverse_gauss_jordan(std::ptrdiff_t n,
                                  const REAL* A_data,
                                  std::ptrdiff_t lda,
                                  std::vector<REAL>& inverse) {
    using A = Rarith<REAL>;

    std::ptrdiff_t width = 2 * n;
    std::vector<REAL> aug(static_cast<std::size_t>(n * width), A::zero());
    inverse.assign(static_cast<std::size_t>(n * n), A::zero());

    for (std::ptrdiff_t row = 0; row < n; ++row) {
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            aug[static_cast<std::size_t>(row * width + col)] = A_data[row * lda + col];
        }
        aug[static_cast<std::size_t>(row * width + n + row)] = A::one();
    }

    for (std::ptrdiff_t col = 0; col < n; ++col) {
        std::ptrdiff_t pivot_row = col;
        REAL pivot_abs = A::abs(aug[static_cast<std::size_t>(col * width + col)]);
        for (std::ptrdiff_t row = col + 1; row < n; ++row) {
            REAL candidate_abs = A::abs(aug[static_cast<std::size_t>(row * width + col)]);
            if (candidate_abs > pivot_abs) {
                pivot_abs = candidate_abs;
                pivot_row = row;
            }
        }
        if (!A::is_finite(pivot_abs) || pivot_abs == A::zero()) {
            return false;
        }
        if (pivot_row != col) {
            for (std::ptrdiff_t j = 0; j < width; ++j) {
                REAL tmp = aug[static_cast<std::size_t>(col * width + j)];
                aug[static_cast<std::size_t>(col * width + j)] = aug[static_cast<std::size_t>(pivot_row * width + j)];
                aug[static_cast<std::size_t>(pivot_row * width + j)] = tmp;
            }
        }

        REAL pivot = aug[static_cast<std::size_t>(col * width + col)];
        if (!A::is_finite(pivot) || pivot == A::zero()) {
            return false;
        }
        for (std::ptrdiff_t j = 0; j < width; ++j) {
            REAL normalized = aug[static_cast<std::size_t>(col * width + j)] / pivot;
            if (!A::is_finite(normalized)) {
                return false;
            }
            aug[static_cast<std::size_t>(col * width + j)] = normalized;
        }

        for (std::ptrdiff_t row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }
            REAL factor = aug[static_cast<std::size_t>(row * width + col)];
            if (!A::is_finite(factor)) {
                return false;
            }
            if (factor == A::zero()) {
                continue;
            }
            for (std::ptrdiff_t j = 0; j < width; ++j) {
                REAL scaled = factor * aug[static_cast<std::size_t>(col * width + j)];
                REAL updated = aug[static_cast<std::size_t>(row * width + j)] - scaled;
                if (!A::is_finite(updated)) {
                    return false;
                }
                aug[static_cast<std::size_t>(row * width + j)] = updated;
            }
        }
    }

    for (std::ptrdiff_t row = 0; row < n; ++row) {
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            inverse[static_cast<std::size_t>(row * n + col)] = aug[static_cast<std::size_t>(row * width + n + col)];
        }
    }
    return finite_contiguous_matrix(n, n, inverse);
}

template <class REAL>
bool compute_solution_from_inverse(std::ptrdiff_t n,
                                   std::ptrdiff_t nrhs,
                                   const std::vector<REAL>& inverse,
                                   const REAL* B_data,
                                   std::ptrdiff_t ldb,
                                   std::vector<REAL>& solution) {
    using A = Rarith<REAL>;

    solution.assign(static_cast<std::size_t>(n * nrhs), A::zero());
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        for (std::ptrdiff_t rhs = 0; rhs < nrhs; ++rhs) {
            REAL acc = A::zero();
            for (std::ptrdiff_t col = 0; col < n; ++col) {
                REAL prod = inverse[static_cast<std::size_t>(row * n + col)] * B_data[col * ldb + rhs];
                REAL next = acc + prod;
                acc = next;
            }
            if (!A::is_finite(acc)) {
                return false;
            }
            solution[static_cast<std::size_t>(row * nrhs + rhs)] = acc;
        }
    }
    return true;
}

template <class REAL>
bool bound_alpha_inverse_residual(std::ptrdiff_t n,
                                  const std::vector<REAL>& inverse,
                                  const REAL* A_data,
                                  std::ptrdiff_t lda,
                                  REAL& alpha) {
    using A = Rarith<REAL>;

    alpha = A::zero();
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        REAL row_sum = A::zero();
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            Rmidrad<REAL> dot = vRdot_apriori_with_reference_fallback(n,
                                                                      inverse.data() + row * n,
                                                                      1,
                                                                      A_data + col,
                                                                      lda);
            if (dot.status != Rstatus::ok || dot.rad < A::zero()) {
                return false;
            }
            REAL delta = (row == col) ? A::one() : A::zero();
            {
                typename A::round_up scope;
                REAL diff = A::zero();
                if (dot.mid >= delta) {
                    diff = dot.mid - delta;
                } else {
                    diff = delta - dot.mid;
                }
                REAL term = diff + dot.rad;
                REAL next = row_sum + term;
                row_sum = next;
            }
            if (!A::is_finite(row_sum)) {
                return false;
            }
        }
        if (row_sum > alpha) {
            alpha = row_sum;
        }
    }
    return A::is_finite(alpha);
}

template <class REAL>
Rmidrad<REAL> residual_component_from_dot(std::ptrdiff_t n,
                                         const REAL* A_row,
                                         const REAL* x_column,
                                         std::ptrdiff_t incx,
                                         REAL rhs) {
    using A = Rarith<REAL>;

    Rmidrad<REAL> dot = vRdot_apriori_with_reference_fallback(n, A_row, 1, x_column, incx);
    if (dot.status != Rstatus::ok || dot.rad < A::zero()) {
        return status_midrad<REAL>(Rstatus::unbounded);
    }

    REAL dot_lo = A::zero();
    {
        typename A::round_down scope;
        dot_lo = dot.mid - dot.rad;
    }
    REAL dot_hi = A::zero();
    {
        typename A::round_up scope;
        dot_hi = dot.mid + dot.rad;
    }

    REAL lo = A::zero();
    {
        typename A::round_down scope;
        lo = rhs - dot_hi;
    }
    REAL hi = A::zero();
    {
        typename A::round_up scope;
        hi = rhs - dot_lo;
    }
    REAL mid = rhs - dot.mid;
    return Rmake_midrad(lo, hi, mid);
}

template <class REAL>
bool bound_beta_from_residuals(std::ptrdiff_t n,
                               const std::vector<REAL>& inverse,
                               const std::vector<Rmidrad<REAL>>& residuals,
                               REAL& beta) {
    using A = Rarith<REAL>;

    std::vector<REAL> residual_abs(static_cast<std::size_t>(n), A::zero());
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        const Rmidrad<REAL>& box = residuals[static_cast<std::size_t>(row)];
        if (box.status != Rstatus::ok || box.rad < A::zero()) {
            return false;
        }
        {
            typename A::round_up scope;
            REAL abs_mid = A::abs(box.mid);
            REAL bound = abs_mid + box.rad;
            residual_abs[static_cast<std::size_t>(row)] = bound;
        }
        if (!A::is_finite(residual_abs[static_cast<std::size_t>(row)])) {
            return false;
        }
    }

    beta = A::zero();
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        REAL row_sum = A::zero();
        {
            typename A::round_up scope;
            for (std::ptrdiff_t col = 0; col < n; ++col) {
                REAL coeff_abs = A::abs(inverse[static_cast<std::size_t>(row * n + col)]);
                REAL prod = coeff_abs * residual_abs[static_cast<std::size_t>(col)];
                REAL next = row_sum + prod;
                row_sum = next;
            }
        }
        if (!A::is_finite(row_sum)) {
            return false;
        }
        if (row_sum > beta) {
            beta = row_sum;
        }
    }
    return A::is_finite(beta);
}

template <class REAL>
bool solution_radius_from_alpha_beta(REAL alpha, REAL beta, REAL& radius) {
    using A = Rarith<REAL>;

    REAL denominator = A::zero();
    {
        typename A::round_down scope;
        denominator = A::one() - alpha;
    }
    if (!A::is_finite(denominator) || denominator <= A::zero()) {
        return false;
    }
    {
        typename A::round_up scope;
        radius = beta / denominator;
    }
    return A::is_finite(radius) && radius >= A::zero();
}

template <class REAL>
bool certify_rhs_column(std::ptrdiff_t n,
                        std::ptrdiff_t nrhs,
                        std::ptrdiff_t rhs,
                        const REAL* A_data,
                        std::ptrdiff_t lda,
                        const REAL* B_data,
                        std::ptrdiff_t ldb,
                        const std::vector<REAL>& inverse,
                        const std::vector<REAL>& solution,
                        REAL alpha,
                        REAL& radius) {
    using A = Rarith<REAL>;

    std::vector<Rmidrad<REAL>> residuals(static_cast<std::size_t>(n));
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        residuals[static_cast<std::size_t>(row)] = residual_component_from_dot(n,
                                                                              A_data + row * lda,
                                                                              solution.data() + rhs,
                                                                              nrhs,
                                                                              B_data[row * ldb + rhs]);
        if (residuals[static_cast<std::size_t>(row)].status != Rstatus::ok) {
            return false;
        }
    }

    REAL beta = A::zero();
    if (!bound_beta_from_residuals(n, inverse, residuals, beta)) {
        return false;
    }
    if (!solution_radius_from_alpha_beta(alpha, beta, radius)) {
        return false;
    }
    return true;
}

template <class REAL>
VerificationStatus vRgesv_matrix_core(std::ptrdiff_t n,
                                      std::ptrdiff_t nrhs,
                                      const REAL* A_data,
                                      std::ptrdiff_t lda,
                                      const REAL* B_data,
                                      std::ptrdiff_t ldb,
                                      Rmidrad<REAL>* X_data,
                                      std::ptrdiff_t ldx) {
    using A = Rarith<REAL>;

    if (!finite_square_matrix(n, A_data, lda) || !finite_rhs_matrix(n, nrhs, B_data, ldb)) {
        fill_matrix_boxes(n, nrhs, X_data, ldx, Rstatus::non_finite);
        return VerificationStatus::Unverified;
    }

    std::vector<REAL> inverse;
    if (!compute_inverse_gauss_jordan(n, A_data, lda, inverse)) {
        fill_matrix_boxes(n, nrhs, X_data, ldx, Rstatus::unbounded);
        return VerificationStatus::Unverified;
    }

    std::vector<REAL> solution;
    if (!compute_solution_from_inverse(n, nrhs, inverse, B_data, ldb, solution)) {
        fill_matrix_boxes(n, nrhs, X_data, ldx, Rstatus::unbounded);
        return VerificationStatus::Unverified;
    }

    REAL alpha = A::zero();
    if (!bound_alpha_inverse_residual(n, inverse, A_data, lda, alpha) || !(alpha < A::one())) {
        fill_matrix_boxes(n, nrhs, X_data, ldx, Rstatus::unbounded);
        return VerificationStatus::Unverified;
    }

    for (std::ptrdiff_t rhs = 0; rhs < nrhs; ++rhs) {
        REAL radius = A::zero();
        if (!certify_rhs_column(n, nrhs, rhs, A_data, lda, B_data, ldb, inverse, solution, alpha, radius)) {
            fill_matrix_boxes(n, nrhs, X_data, ldx, Rstatus::unbounded);
            return VerificationStatus::Unverified;
        }
        for (std::ptrdiff_t row = 0; row < n; ++row) {
            REAL mid = solution[static_cast<std::size_t>(row * nrhs + rhs)];
            if (!A::is_finite(mid)) {
                fill_matrix_boxes(n, nrhs, X_data, ldx, Rstatus::unbounded);
                return VerificationStatus::Unverified;
            }
            X_data[row * ldx + rhs] = {mid, radius, Rstatus::ok};
        }
    }

    return VerificationStatus::Verified;
}

} // namespace detail

template <class REAL>
Rmidrad<REAL> vRsum(std::ptrdiff_t n, const REAL* x, std::ptrdiff_t incx) {
    using A = Rarith<REAL>;

    if (n < 0) {
        return detail::invalid_midrad<REAL>();
    }
    if (n == 0) {
        return {A::zero(), A::zero(), Rstatus::ok};
    }
    if (x == nullptr || incx < 1) {
        return detail::invalid_midrad<REAL>();
    }

    for (std::ptrdiff_t i = 0; i < n; ++i) {
        if (!A::is_finite(x[i * incx])) {
            return detail::non_finite_midrad<REAL>();
        }
    }

    REAL sup = A::zero();
    {
        typename A::round_up scope;
        REAL acc = A::zero();
        for (std::ptrdiff_t i = 0; i < n; ++i) {
            REAL next = acc + x[i * incx];
            acc = next;
        }
        sup = acc;
    }

    REAL inf = A::zero();
    {
        typename A::round_down scope;
        REAL acc = A::zero();
        for (std::ptrdiff_t i = 0; i < n; ++i) {
            REAL next = acc + x[i * incx];
            acc = next;
        }
        inf = acc;
    }

    REAL mid = Rsum(n, x, incx);
    return Rmake_midrad(inf, sup, mid);
}

template <class REAL>
Rmidrad<REAL> vRdot(std::ptrdiff_t n,
                    const REAL* x,
                    std::ptrdiff_t incx,
                    const REAL* y,
                    std::ptrdiff_t incy) {
    using A = Rarith<REAL>;

    if (n < 0) {
        return detail::invalid_midrad<REAL>();
    }
    if (n == 0) {
        return {A::zero(), A::zero(), Rstatus::ok};
    }
    if (x == nullptr || y == nullptr || incx < 1 || incy < 1) {
        return detail::invalid_midrad<REAL>();
    }

    for (std::ptrdiff_t i = 0; i < n; ++i) {
        if (!A::is_finite(x[i * incx]) || !A::is_finite(y[i * incy])) {
            return detail::non_finite_midrad<REAL>();
        }
    }

    REAL sup = A::zero();
    {
        typename A::round_up scope;
        REAL acc = A::zero();
        for (std::ptrdiff_t i = 0; i < n; ++i) {
            REAL prod = x[i * incx] * y[i * incy];
            REAL next = acc + prod;
            acc = next;
        }
        sup = acc;
    }

    REAL inf = A::zero();
    {
        typename A::round_down scope;
        REAL acc = A::zero();
        for (std::ptrdiff_t i = 0; i < n; ++i) {
            REAL prod = x[i * incx] * y[i * incy];
            REAL next = acc + prod;
            acc = next;
        }
        inf = acc;
    }

    REAL mid = Rdot(n, x, incx, y, incy);
    return Rmake_midrad(inf, sup, mid);
}


template <class REAL>
Rmidrad<REAL> vRdot_apriori(std::ptrdiff_t n,
                            const REAL* x,
                            std::ptrdiff_t incx,
                            const REAL* y,
                            std::ptrdiff_t incy) {
    using A = Rarith<REAL>;

    if (n < 0) {
        return detail::invalid_midrad<REAL>();
    }
    if (n == 0) {
        return {A::zero(), A::zero(), Rstatus::ok};
    }
    if (x == nullptr || y == nullptr || incx < 1 || incy < 1) {
        return detail::invalid_midrad<REAL>();
    }

    for (std::ptrdiff_t i = 0; i < n; ++i) {
        if (!A::is_finite(x[i * incx]) || !A::is_finite(y[i * incy])) {
            return detail::non_finite_midrad<REAL>();
        }
    }

    detail::clear_underflow_flag<REAL>();
    REAL mid = Rdot(n, x, incx, y, incy);
    if (detail::underflow_flag_raised<REAL>() || !A::is_finite(mid)) {
        return detail::unbounded_midrad_with_mid(mid);
    }

    if (!detail::apriori_length_gate(n, A::precision_bits())) {
        return detail::unbounded_midrad_with_mid(mid);
    }

    detail::clear_underflow_flag<REAL>();
    REAL S_up = Rsum_abs_dot_upper(n, x, incx, y, incy);
    if (detail::underflow_flag_raised<REAL>() || !A::is_finite(S_up)) {
        return detail::unbounded_midrad_with_mid(mid);
    }

    detail::clear_underflow_flag<REAL>();
    REAL rad = detail::Rdot_apriori_radius(n, mid, S_up);
    if (detail::underflow_flag_raised<REAL>() || !A::is_finite(rad)) {
        return detail::unbounded_midrad_with_mid(mid);
    }

    return {mid, rad, Rstatus::ok};
}


template <class REAL>
Rstatus vRresidual(std::ptrdiff_t m,
                   std::ptrdiff_t n,
                   const REAL* A_data,
                   std::ptrdiff_t lda,
                   const REAL* x,
                   const REAL* b,
                   Rmidrad<REAL>* out) {
    using Arith = Rarith<REAL>;

    if (m < 0 || n < 0) {
        return Rstatus::invalid_input;
    }
    if (m == 0) {
        return Rstatus::ok;
    }
    if (out == nullptr || b == nullptr) {
        return Rstatus::invalid_input;
    }
    if (n == 0) {
        Rstatus worst = Rstatus::ok;
        for (std::ptrdiff_t i = 0; i < m; ++i) {
            if (!Arith::is_finite(b[i])) {
                out[i] = detail::non_finite_midrad<REAL>();
            } else {
                out[i] = {b[i], Arith::zero(), Rstatus::ok};
            }
            worst = detail::worst_status(worst, out[i].status);
        }
        return worst;
    }
    if (A_data == nullptr || x == nullptr || lda < n) {
        return Rstatus::invalid_input;
    }

    Rstatus worst = Rstatus::ok;
    for (std::ptrdiff_t row = 0; row < m; ++row) {
        bool finite = Arith::is_finite(b[row]);
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            if (!Arith::is_finite(A_data[row * lda + col]) || !Arith::is_finite(x[col])) {
                finite = false;
            }
        }
        if (!finite) {
            out[row] = detail::non_finite_midrad<REAL>();
            worst = detail::worst_status(worst, out[row].status);
            continue;
        }

        REAL sup = Arith::zero();
        {
            typename Arith::round_up scope;
            REAL acc = b[row];
            for (std::ptrdiff_t col = 0; col < n; ++col) {
                REAL na = -A_data[row * lda + col];
                REAL prod = na * x[col];
                REAL next = acc + prod;
                acc = next;
            }
            sup = acc;
        }

        REAL inf = Arith::zero();
        {
            typename Arith::round_down scope;
            REAL acc = b[row];
            for (std::ptrdiff_t col = 0; col < n; ++col) {
                REAL na = -A_data[row * lda + col];
                REAL prod = na * x[col];
                REAL next = acc + prod;
                acc = next;
            }
            inf = acc;
        }

        REAL mid = Rmidpoint(inf, sup);
        out[row] = Rmake_midrad(inf, sup, mid);
        worst = detail::worst_status(worst, out[row].status);
    }

    return worst;
}


template <class REAL>
VerificationStatus vRgemv_point(std::ptrdiff_t m,
                                std::ptrdiff_t n,
                                const REAL* A_data,
                                std::ptrdiff_t lda,
                                const REAL* x,
                                std::ptrdiff_t incx,
                                Rmidrad<REAL>* out) {
    using Arith = Rarith<REAL>;

    if (m < 0 || n < 0) {
        return VerificationStatus::InvalidInput;
    }
    if (m == 0) {
        return VerificationStatus::Verified;
    }
    if (out == nullptr) {
        return VerificationStatus::InvalidInput;
    }
    if (n == 0) {
        for (std::ptrdiff_t row = 0; row < m; ++row) {
            out[row] = {Arith::zero(), Arith::zero(), Rstatus::ok};
        }
        return VerificationStatus::Verified;
    }
    if (A_data == nullptr || x == nullptr || lda < n || incx < 1) {
        return VerificationStatus::InvalidInput;
    }

    VerificationStatus status = VerificationStatus::Verified;
    for (std::ptrdiff_t row = 0; row < m; ++row) {
        Rmidrad<REAL> box = detail::vRdot_apriori_with_reference_fallback(n, A_data + row * lda, 1, x, incx);
        out[row] = box;
        status = detail::worst_verification_status(status, detail::verification_status_from_box(box.status));
    }
    return status;
}


template <class REAL>
VerificationStatus vRgemm_point(std::ptrdiff_t m,
                                std::ptrdiff_t n,
                                std::ptrdiff_t k,
                                const REAL* A_data,
                                std::ptrdiff_t lda,
                                const REAL* B_data,
                                std::ptrdiff_t ldb,
                                Rmidrad<REAL>* C_data,
                                std::ptrdiff_t ldc) {
    using Arith = Rarith<REAL>;

    if (m < 0 || n < 0 || k < 0) {
        return VerificationStatus::InvalidInput;
    }
    if (m == 0 || n == 0) {
        return VerificationStatus::Verified;
    }
    if (C_data == nullptr || ldc < n) {
        return VerificationStatus::InvalidInput;
    }
    if (k == 0) {
        for (std::ptrdiff_t row = 0; row < m; ++row) {
            for (std::ptrdiff_t col = 0; col < n; ++col) {
                C_data[row * ldc + col] = {Arith::zero(), Arith::zero(), Rstatus::ok};
            }
        }
        return VerificationStatus::Verified;
    }
    if (A_data == nullptr || B_data == nullptr || lda < k || ldb < n) {
        return VerificationStatus::InvalidInput;
    }

    VerificationStatus status = VerificationStatus::Verified;
    for (std::ptrdiff_t row = 0; row < m; ++row) {
        for (std::ptrdiff_t col = 0; col < n; ++col) {
            Rmidrad<REAL> box = detail::vRdot_apriori_with_reference_fallback(k, A_data + row * lda, 1, B_data + col, ldb);
            C_data[row * ldc + col] = box;
            status = detail::worst_verification_status(status, detail::verification_status_from_box(box.status));
        }
    }
    return status;
}


template <class REAL>
VerificationStatus vRgesv(std::ptrdiff_t n,
                          const REAL* A_data,
                          std::ptrdiff_t lda,
                          const REAL* b,
                          std::ptrdiff_t incb,
                          Rmidrad<REAL>* x,
                          std::ptrdiff_t incx) {
    using A = Rarith<REAL>;

    if (n < 0) {
        return VerificationStatus::InvalidInput;
    }
    if (n == 0) {
        return VerificationStatus::Verified;
    }
    if (A_data == nullptr || b == nullptr || x == nullptr || lda < n || incb < 1 || incx < 1) {
        return VerificationStatus::InvalidInput;
    }

    std::ptrdiff_t A_last = detail::matrix_last_index(n, n, lda);
    std::ptrdiff_t b_last = detail::strided_last_index(n, incb);
    std::ptrdiff_t x_last = detail::strided_last_index(n, incx);
    if (detail::storage_ranges_overlap(A_data, A_last, x, x_last) ||
        detail::storage_ranges_overlap(b, b_last, x, x_last)) {
        return VerificationStatus::InvalidInput;
    }

    std::vector<REAL> compact_b(static_cast<std::size_t>(n), A::zero());
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        compact_b[static_cast<std::size_t>(row)] = b[row * incb];
    }
    std::vector<Rmidrad<REAL>> compact_x(static_cast<std::size_t>(n));
    VerificationStatus status = vRgesv(n, 1, A_data, lda, compact_b.data(), 1, compact_x.data(), 1);
    for (std::ptrdiff_t row = 0; row < n; ++row) {
        x[row * incx] = compact_x[static_cast<std::size_t>(row)];
    }
    return status;
}

template <class REAL>
VerificationStatus vRgesv(std::ptrdiff_t n,
                          std::ptrdiff_t nrhs,
                          const REAL* A_data,
                          std::ptrdiff_t lda,
                          const REAL* B_data,
                          std::ptrdiff_t ldb,
                          Rmidrad<REAL>* X_data,
                          std::ptrdiff_t ldx) {
    if (n < 0 || nrhs < 0) {
        return VerificationStatus::InvalidInput;
    }
    if (n == 0 || nrhs == 0) {
        return VerificationStatus::Verified;
    }
    if (A_data == nullptr || B_data == nullptr || X_data == nullptr || lda < n || ldb < nrhs || ldx < nrhs) {
        return VerificationStatus::InvalidInput;
    }

    std::ptrdiff_t A_last = detail::matrix_last_index(n, n, lda);
    std::ptrdiff_t B_last = detail::matrix_last_index(n, nrhs, ldb);
    std::ptrdiff_t X_last = detail::matrix_last_index(n, nrhs, ldx);
    if (detail::storage_ranges_overlap(A_data, A_last, X_data, X_last) ||
        detail::storage_ranges_overlap(B_data, B_last, X_data, X_last)) {
        return VerificationStatus::InvalidInput;
    }

    return detail::vRgesv_matrix_core(n, nrhs, A_data, lda, B_data, ldb, X_data, ldx);
}

template <class REAL>
VerificationStatus vRgeinv(std::ptrdiff_t n,
                           const REAL* A_data,
                           std::ptrdiff_t lda,
                           Rmidrad<REAL>* Ainv_data,
                           std::ptrdiff_t ldinv) {
    using A = Rarith<REAL>;

    if (n < 0) {
        return VerificationStatus::InvalidInput;
    }
    if (n == 0) {
        return VerificationStatus::Verified;
    }
    if (A_data == nullptr || Ainv_data == nullptr || lda < n || ldinv < n) {
        return VerificationStatus::InvalidInput;
    }

    std::ptrdiff_t A_last = detail::matrix_last_index(n, n, lda);
    std::ptrdiff_t Ainv_last = detail::matrix_last_index(n, n, ldinv);
    if (detail::storage_ranges_overlap(A_data, A_last, Ainv_data, Ainv_last)) {
        return VerificationStatus::InvalidInput;
    }

    std::vector<REAL> identity(static_cast<std::size_t>(n * n), A::zero());
    for (std::ptrdiff_t i = 0; i < n; ++i) {
        identity[static_cast<std::size_t>(i * n + i)] = A::one();
    }

    return vRgesv(n, n, A_data, lda, identity.data(), n, Ainv_data, ldinv);
}

template <class REAL>
Rcond_bounds<REAL> vRgecon(std::ptrdiff_t n, const REAL* A_data, std::ptrdiff_t lda) {
    using A = Rarith<REAL>;

    if (n < 0) {
        return detail::invalid_cond_bounds<REAL>();
    }
    if (n == 0) {
        return detail::empty_cond_bounds<REAL>();
    }
    if (A_data == nullptr || lda < n) {
        return detail::invalid_cond_bounds<REAL>();
    }
    if (!detail::finite_square_matrix(n, A_data, lda)) {
        return detail::non_finite_cond_bounds<REAL>();
    }

    REAL normA_upper = A::zero();
    if (!detail::infinity_norm_point_upper(n, A_data, lda, normA_upper)) {
        return detail::unbounded_cond_bounds<REAL>(normA_upper);
    }

    std::vector<Rmidrad<REAL>> Ainv_data(static_cast<std::size_t>(n * n));
    VerificationStatus inverse_status = vRgeinv(n, A_data, lda, Ainv_data.data(), n);
    if (inverse_status == VerificationStatus::InvalidInput) {
        return detail::invalid_cond_bounds<REAL>();
    }
    if (inverse_status != VerificationStatus::Verified) {
        return detail::unbounded_cond_bounds<REAL>(normA_upper);
    }

    Rstatus inverse_box_status = detail::worst_midrad_matrix_status(n, Ainv_data.data(), n);
    if (inverse_box_status != Rstatus::ok) {
        return detail::unbounded_cond_bounds<REAL>(normA_upper);
    }

    REAL normAinv_upper = A::zero();
    if (!detail::infinity_norm_midrad_upper(n, Ainv_data.data(), n, normAinv_upper)) {
        return detail::unbounded_cond_bounds<REAL>(normA_upper);
    }

    REAL kappa_upper = A::zero();
    {
        typename A::round_up scope;
        REAL product = normA_upper * normAinv_upper;
        kappa_upper = product;
    }
    if (!A::is_finite(kappa_upper) || !(kappa_upper > A::zero())) {
        return detail::unbounded_cond_bounds<REAL>(normA_upper);
    }

    REAL rcond_lower = A::zero();
    {
        typename A::round_down scope;
        REAL reciprocal = A::one() / kappa_upper;
        rcond_lower = reciprocal;
    }
    if (!A::is_finite(rcond_lower) || rcond_lower < A::zero()) {
        return detail::unbounded_cond_bounds<REAL>(normA_upper);
    }

    return {normA_upper, normAinv_upper, kappa_upper, rcond_lower, Rstatus::ok};
}

extern template float Rsum<float>(std::ptrdiff_t, const float*, std::ptrdiff_t);
extern template double Rsum<double>(std::ptrdiff_t, const double*, std::ptrdiff_t);
extern template float Rdot<float>(std::ptrdiff_t, const float*, std::ptrdiff_t, const float*, std::ptrdiff_t);
extern template double Rdot<double>(std::ptrdiff_t, const double*, std::ptrdiff_t, const double*, std::ptrdiff_t);
extern template float Rsum_abs_dot_upper<float>(std::ptrdiff_t,
                                               const float*,
                                               std::ptrdiff_t,
                                               const float*,
                                               std::ptrdiff_t);
extern template double Rsum_abs_dot_upper<double>(std::ptrdiff_t,
                                                 const double*,
                                                 std::ptrdiff_t,
                                                 const double*,
                                                 std::ptrdiff_t);
extern template Rmidrad<float> vRsum<float>(std::ptrdiff_t, const float*, std::ptrdiff_t);
extern template Rmidrad<double> vRsum<double>(std::ptrdiff_t, const double*, std::ptrdiff_t);
extern template Rmidrad<float> vRdot<float>(std::ptrdiff_t, const float*, std::ptrdiff_t, const float*, std::ptrdiff_t);
extern template Rmidrad<double> vRdot<double>(std::ptrdiff_t, const double*, std::ptrdiff_t, const double*, std::ptrdiff_t);
extern template Rmidrad<float> vRdot_apriori<float>(std::ptrdiff_t,
                                                   const float*,
                                                   std::ptrdiff_t,
                                                   const float*,
                                                   std::ptrdiff_t);
extern template Rmidrad<double> vRdot_apriori<double>(std::ptrdiff_t,
                                                     const double*,
                                                     std::ptrdiff_t,
                                                     const double*,
                                                     std::ptrdiff_t);
extern template Rstatus vRresidual<float>(std::ptrdiff_t,
                                          std::ptrdiff_t,
                                          const float*,
                                          std::ptrdiff_t,
                                          const float*,
                                          const float*,
                                          Rmidrad<float>*);
extern template Rstatus vRresidual<double>(std::ptrdiff_t,
                                           std::ptrdiff_t,
                                           const double*,
                                           std::ptrdiff_t,
                                           const double*,
                                           const double*,
                                           Rmidrad<double>*);
extern template VerificationStatus vRgemv_point<float>(std::ptrdiff_t,
                                                       std::ptrdiff_t,
                                                       const float*,
                                                       std::ptrdiff_t,
                                                       const float*,
                                                       std::ptrdiff_t,
                                                       Rmidrad<float>*);
extern template VerificationStatus vRgemv_point<double>(std::ptrdiff_t,
                                                        std::ptrdiff_t,
                                                        const double*,
                                                        std::ptrdiff_t,
                                                        const double*,
                                                        std::ptrdiff_t,
                                                        Rmidrad<double>*);
extern template VerificationStatus vRgemm_point<float>(std::ptrdiff_t,
                                                       std::ptrdiff_t,
                                                       std::ptrdiff_t,
                                                       const float*,
                                                       std::ptrdiff_t,
                                                       const float*,
                                                       std::ptrdiff_t,
                                                       Rmidrad<float>*,
                                                       std::ptrdiff_t);
extern template VerificationStatus vRgemm_point<double>(std::ptrdiff_t,
                                                        std::ptrdiff_t,
                                                        std::ptrdiff_t,
                                                        const double*,
                                                        std::ptrdiff_t,
                                                        const double*,
                                                        std::ptrdiff_t,
                                                        Rmidrad<double>*,
                                                        std::ptrdiff_t);
extern template VerificationStatus vRgesv<float>(std::ptrdiff_t,
                                                const float*,
                                                std::ptrdiff_t,
                                                const float*,
                                                std::ptrdiff_t,
                                                Rmidrad<float>*,
                                                std::ptrdiff_t);
extern template VerificationStatus vRgesv<double>(std::ptrdiff_t,
                                                 const double*,
                                                 std::ptrdiff_t,
                                                 const double*,
                                                 std::ptrdiff_t,
                                                 Rmidrad<double>*,
                                                 std::ptrdiff_t);
extern template VerificationStatus vRgesv<float>(std::ptrdiff_t,
                                                std::ptrdiff_t,
                                                const float*,
                                                std::ptrdiff_t,
                                                const float*,
                                                std::ptrdiff_t,
                                                Rmidrad<float>*,
                                                std::ptrdiff_t);
extern template VerificationStatus vRgesv<double>(std::ptrdiff_t,
                                                 std::ptrdiff_t,
                                                 const double*,
                                                 std::ptrdiff_t,
                                                 const double*,
                                                 std::ptrdiff_t,
                                                 Rmidrad<double>*,
                                                 std::ptrdiff_t);
extern template VerificationStatus vRgeinv<float>(std::ptrdiff_t,
                                                  const float*,
                                                  std::ptrdiff_t,
                                                  Rmidrad<float>*,
                                                  std::ptrdiff_t);
extern template VerificationStatus vRgeinv<double>(std::ptrdiff_t,
                                                   const double*,
                                                   std::ptrdiff_t,
                                                   Rmidrad<double>*,
                                                   std::ptrdiff_t);

extern template Rcond_bounds<float> vRgecon<float>(std::ptrdiff_t, const float*, std::ptrdiff_t);
extern template Rcond_bounds<double> vRgecon<double>(std::ptrdiff_t, const double*, std::ptrdiff_t);

#ifdef VMPLAPACK_ENABLE_MPFR
extern template mpfrxx::mpfr_class Rsum<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                            const mpfrxx::mpfr_class*,
                                                            std::ptrdiff_t);
extern template mpfrxx::mpfr_class Rdot<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                            const mpfrxx::mpfr_class*,
                                                            std::ptrdiff_t,
                                                            const mpfrxx::mpfr_class*,
                                                            std::ptrdiff_t);
extern template mpfrxx::mpfr_class Rsum_abs_dot_upper<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                                          const mpfrxx::mpfr_class*,
                                                                          std::ptrdiff_t,
                                                                          const mpfrxx::mpfr_class*,
                                                                          std::ptrdiff_t);
extern template Rmidrad<mpfrxx::mpfr_class> vRsum<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                                      const mpfrxx::mpfr_class*,
                                                                      std::ptrdiff_t);
extern template Rmidrad<mpfrxx::mpfr_class> vRdot<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                                      const mpfrxx::mpfr_class*,
                                                                      std::ptrdiff_t,
                                                                      const mpfrxx::mpfr_class*,
                                                                      std::ptrdiff_t);
extern template Rmidrad<mpfrxx::mpfr_class> vRdot_apriori<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                                              const mpfrxx::mpfr_class*,
                                                                              std::ptrdiff_t,
                                                                              const mpfrxx::mpfr_class*,
                                                                              std::ptrdiff_t);
extern template Rstatus vRresidual<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                       std::ptrdiff_t,
                                                       const mpfrxx::mpfr_class*,
                                                       std::ptrdiff_t,
                                                       const mpfrxx::mpfr_class*,
                                                       const mpfrxx::mpfr_class*,
                                                       Rmidrad<mpfrxx::mpfr_class>*);
extern template VerificationStatus vRgemv_point<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                                    std::ptrdiff_t,
                                                                    const mpfrxx::mpfr_class*,
                                                                    std::ptrdiff_t,
                                                                    const mpfrxx::mpfr_class*,
                                                                    std::ptrdiff_t,
                                                                    Rmidrad<mpfrxx::mpfr_class>*);
extern template VerificationStatus vRgemm_point<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                                    std::ptrdiff_t,
                                                                    std::ptrdiff_t,
                                                                    const mpfrxx::mpfr_class*,
                                                                    std::ptrdiff_t,
                                                                    const mpfrxx::mpfr_class*,
                                                                    std::ptrdiff_t,
                                                                    Rmidrad<mpfrxx::mpfr_class>*,
                                                                    std::ptrdiff_t);
extern template VerificationStatus vRgesv<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                             const mpfrxx::mpfr_class*,
                                                             std::ptrdiff_t,
                                                             const mpfrxx::mpfr_class*,
                                                             std::ptrdiff_t,
                                                             Rmidrad<mpfrxx::mpfr_class>*,
                                                             std::ptrdiff_t);
extern template VerificationStatus vRgesv<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                             std::ptrdiff_t,
                                                             const mpfrxx::mpfr_class*,
                                                             std::ptrdiff_t,
                                                             const mpfrxx::mpfr_class*,
                                                             std::ptrdiff_t,
                                                             Rmidrad<mpfrxx::mpfr_class>*,
                                                             std::ptrdiff_t);
extern template VerificationStatus vRgeinv<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                               const mpfrxx::mpfr_class*,
                                                               std::ptrdiff_t,
                                                               Rmidrad<mpfrxx::mpfr_class>*,
                                                               std::ptrdiff_t);
extern template Rcond_bounds<mpfrxx::mpfr_class> vRgecon<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                                            const mpfrxx::mpfr_class*,
                                                                            std::ptrdiff_t);
#endif

} // namespace vmplapack
