// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <cstddef>

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
Rmidrad<REAL> vRsum(std::ptrdiff_t n, const REAL* x, std::ptrdiff_t incx = 1);

template <class REAL>
Rmidrad<REAL> vRdot(std::ptrdiff_t n,
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

} // namespace vmplapack
