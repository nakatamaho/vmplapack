// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

namespace vmplapack {

template VerificationStatus vRgeinv<float>(std::ptrdiff_t,
                                           const float*,
                                           std::ptrdiff_t,
                                           Rmidrad<float>*,
                                           std::ptrdiff_t);
template VerificationStatus vRgeinv<double>(std::ptrdiff_t,
                                            const double*,
                                            std::ptrdiff_t,
                                            Rmidrad<double>*,
                                            std::ptrdiff_t);

#ifdef VMPLAPACK_ENABLE_MPFR
template VerificationStatus vRgeinv<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                        const mpfrxx::mpfr_class*,
                                                        std::ptrdiff_t,
                                                        Rmidrad<mpfrxx::mpfr_class>*,
                                                        std::ptrdiff_t);
#endif

} // namespace vmplapack
