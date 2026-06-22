// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

namespace vmplapack {

template Rstatus vRresidual<float>(std::ptrdiff_t,
                                   std::ptrdiff_t,
                                   const float*,
                                   std::ptrdiff_t,
                                   const float*,
                                   const float*,
                                   Rmidrad<float>*);
template Rstatus vRresidual<double>(std::ptrdiff_t,
                                    std::ptrdiff_t,
                                    const double*,
                                    std::ptrdiff_t,
                                    const double*,
                                    const double*,
                                    Rmidrad<double>*);

#ifdef VMPLAPACK_ENABLE_MPFR
template Rstatus vRresidual<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                std::ptrdiff_t,
                                                const mpfrxx::mpfr_class*,
                                                std::ptrdiff_t,
                                                const mpfrxx::mpfr_class*,
                                                const mpfrxx::mpfr_class*,
                                                Rmidrad<mpfrxx::mpfr_class>*);
#endif

} // namespace vmplapack
