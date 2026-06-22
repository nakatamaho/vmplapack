// Copyright (c) 2026, NAKATA Maho
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <cassert>
#include <cstddef>

#include <vmplapack/vmplapack_arith.h>

namespace vmplapack {

enum class Rstatus {
    ok,
    unbounded,
    non_finite,
    invalid_input
};

inline int Rstatus_rank(Rstatus s) {
    switch (s) {
        case Rstatus::ok:
            return 0;
        case Rstatus::unbounded:
            return 1;
        case Rstatus::non_finite:
            return 2;
        case Rstatus::invalid_input:
            return 3;
    }
    return 3;
}

template <class REAL>
struct Rmidrad {
    REAL mid;
    REAL rad;
    Rstatus status;
};

template <class REAL>
REAL Rmidpoint(REAL lo, REAL hi) {
    REAL h = Rarith<REAL>::half();
    REAL a = lo * h;
    REAL b = hi * h;
    REAL m = a + b;
    return m;
}

template <class REAL>
REAL Rupward_abs_diff(REAL a, REAL b) {
    typename Rarith<REAL>::round_up g;
    REAL d = (a >= b) ? (a - b) : (b - a);
    return d;
}

template <class REAL>
Rmidrad<REAL> Rmake_midrad(REAL lo, REAL hi, REAL mid) {
    using A = Rarith<REAL>;

    if (!A::is_finite(lo) || !A::is_finite(hi)) {
        return {A::zero(), A::infinity(), Rstatus::unbounded};
    }
    if (hi < lo) {
        assert(false && "Rmake_midrad: lo > hi");
        return {A::zero(), A::infinity(), Rstatus::unbounded};
    }

    REAL m = mid;
    if (!A::is_finite(m)) {
        m = Rmidpoint(lo, hi);
    }

    REAL r1 = Rupward_abs_diff(m, lo);
    REAL r2 = Rupward_abs_diff(hi, m);
    REAL rad = (r1 > r2) ? r1 : r2;
    if (!A::is_finite(rad)) {
        return {A::zero(), A::infinity(), Rstatus::unbounded};
    }

    return {m, rad, Rstatus::ok};
}

} // namespace vmplapack
