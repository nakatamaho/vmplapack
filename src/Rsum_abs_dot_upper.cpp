// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

namespace vmplapack {

template float Rsum_abs_dot_upper<float>(std::ptrdiff_t, const float*, std::ptrdiff_t, const float*, std::ptrdiff_t);
template double Rsum_abs_dot_upper<double>(std::ptrdiff_t, const double*, std::ptrdiff_t, const double*, std::ptrdiff_t);

#ifdef VMPLAPACK_ENABLE_MPFR
template mpfrxx::mpfr_class Rsum_abs_dot_upper<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                                  const mpfrxx::mpfr_class*,
                                                                  std::ptrdiff_t,
                                                                  const mpfrxx::mpfr_class*,
                                                                  std::ptrdiff_t);
#endif

} // namespace vmplapack
