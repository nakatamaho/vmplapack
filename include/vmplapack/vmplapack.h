// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#if defined(__FAST_MATH__)
#error "vMPLAPACK requires IEEE-754 semantics; do not compile with -ffast-math."
#endif

#include <cfloat>

static_assert(FLT_EVAL_METHOD == 0,
              "Excess precision breaks EFT, interval bounds, and the float tier.");

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif
#pragma STDC FENV_ACCESS ON
#pragma STDC FP_CONTRACT OFF
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#include <vmplapack/vmplapack_arith.h>
#include <vmplapack/vmplapack_eft.h>
#include <vmplapack/vmplapack_utils.h>
#include <vmplapack/vmpblas.h>
