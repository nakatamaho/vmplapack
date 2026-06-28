// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

namespace vmplapack {

template VerificationStatus vRgesv<float>(std::ptrdiff_t,
                                          const float*,
                                          std::ptrdiff_t,
                                          const float*,
                                          std::ptrdiff_t,
                                          Rmidrad<float>*,
                                          std::ptrdiff_t);
template VerificationStatus vRgesv<double>(std::ptrdiff_t,
                                           const double*,
                                           std::ptrdiff_t,
                                           const double*,
                                           std::ptrdiff_t,
                                           Rmidrad<double>*,
                                           std::ptrdiff_t);
template VerificationStatus vRgesv<float>(std::ptrdiff_t,
                                          std::ptrdiff_t,
                                          const float*,
                                          std::ptrdiff_t,
                                          const float*,
                                          std::ptrdiff_t,
                                          Rmidrad<float>*,
                                          std::ptrdiff_t);
template VerificationStatus vRgesv<double>(std::ptrdiff_t,
                                           std::ptrdiff_t,
                                           const double*,
                                           std::ptrdiff_t,
                                           const double*,
                                           std::ptrdiff_t,
                                           Rmidrad<double>*,
                                           std::ptrdiff_t);

#ifdef VMPLAPACK_ENABLE_MPFR
template VerificationStatus vRgesv<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                       const mpfrxx::mpfr_class*,
                                                       std::ptrdiff_t,
                                                       const mpfrxx::mpfr_class*,
                                                       std::ptrdiff_t,
                                                       Rmidrad<mpfrxx::mpfr_class>*,
                                                       std::ptrdiff_t);
template VerificationStatus vRgesv<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                       std::ptrdiff_t,
                                                       const mpfrxx::mpfr_class*,
                                                       std::ptrdiff_t,
                                                       const mpfrxx::mpfr_class*,
                                                       std::ptrdiff_t,
                                                       Rmidrad<mpfrxx::mpfr_class>*,
                                                       std::ptrdiff_t);
#endif

} // namespace vmplapack
