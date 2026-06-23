// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#include <vmplapack/vmplapack.h>

namespace vmplapack {

template Rmidrad<float> vRdot_apriori<float>(std::ptrdiff_t, const float*, std::ptrdiff_t, const float*, std::ptrdiff_t);
template Rmidrad<double> vRdot_apriori<double>(std::ptrdiff_t, const double*, std::ptrdiff_t, const double*, std::ptrdiff_t);

#ifdef VMPLAPACK_ENABLE_MPFR
template Rmidrad<mpfrxx::mpfr_class> vRdot_apriori<mpfrxx::mpfr_class>(std::ptrdiff_t,
                                                                      const mpfrxx::mpfr_class*,
                                                                      std::ptrdiff_t,
                                                                      const mpfrxx::mpfr_class*,
                                                                      std::ptrdiff_t);
#endif

} // namespace vmplapack
