// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

namespace vmplapack {

template Rcond_bounds<float> vRgecon<float>(std::ptrdiff_t, const float*, std::ptrdiff_t);
template Rcond_bounds<double> vRgecon<double>(std::ptrdiff_t, const double*, std::ptrdiff_t);

#ifdef VMPLAPACK_ENABLE_MPFR
template Rcond_bounds<mpfrxx::mpfr_class> vRgecon<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                                     const mpfrxx::mpfr_class*,
                                                                     std::ptrdiff_t);
#endif

} // namespace vmplapack
