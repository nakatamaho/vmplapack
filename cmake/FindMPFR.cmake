# Copyright (c) 2026, NAKATA Maho
# SPDX-License-Identifier: BSD-2-Clause

include(FindPackageHandleStandardArgs)

set(_MPFR_HINTS)
foreach(_MPFR_ROOT IN ITEMS
        "${MPFR_ROOT}"
        "$ENV{MPFR_ROOT}"
        "/home/docker/gmpfrxx_mkII/build_release/_deps/gmpfrxx_mkii/MPFR"
        "/home/docker/mplapack/external/i/MPFR")
    if(_MPFR_ROOT)
        list(APPEND _MPFR_HINTS "${_MPFR_ROOT}")
    endif()
endforeach()

find_path(MPFR_INCLUDE_DIR
    NAMES mpfr.h
    HINTS ${_MPFR_HINTS}
    PATH_SUFFIXES include)

find_library(MPFR_LIBRARY
    NAMES mpfr
    HINTS ${_MPFR_HINTS}
    PATH_SUFFIXES lib lib64)

find_package_handle_standard_args(MPFR
    REQUIRED_VARS MPFR_LIBRARY MPFR_INCLUDE_DIR)

if(MPFR_FOUND AND NOT TARGET MPFR::MPFR)
    add_library(MPFR::MPFR UNKNOWN IMPORTED)
    set_target_properties(MPFR::MPFR
        PROPERTIES
            IMPORTED_LOCATION "${MPFR_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDE_DIR}")
    if(TARGET GMP::GMP)
        set_target_properties(MPFR::MPFR
            PROPERTIES
                INTERFACE_LINK_LIBRARIES GMP::GMP)
    endif()
endif()

mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARY)
