# Copyright (c) 2026, NAKATA Maho
# SPDX-License-Identifier: BSD-2-Clause

include(FindPackageHandleStandardArgs)

set(_GMP_HINTS)
foreach(_GMP_ROOT IN ITEMS
        "${GMP_ROOT}"
        "$ENV{GMP_ROOT}"
        "/home/docker/gmpfrxx_mkII/build_release/_deps/gmpfrxx_mkii/GMP"
        "/home/docker/mplapack/external/i/GMP")
    if(_GMP_ROOT)
        list(APPEND _GMP_HINTS "${_GMP_ROOT}")
    endif()
endforeach()

find_path(GMP_INCLUDE_DIR
    NAMES gmp.h
    HINTS ${_GMP_HINTS}
    PATH_SUFFIXES include)

find_library(GMP_LIBRARY
    NAMES gmp
    HINTS ${_GMP_HINTS}
    PATH_SUFFIXES lib lib64)

find_package_handle_standard_args(GMP
    REQUIRED_VARS GMP_LIBRARY GMP_INCLUDE_DIR)

if(GMP_FOUND AND NOT TARGET GMP::GMP)
    add_library(GMP::GMP UNKNOWN IMPORTED)
    set_target_properties(GMP::GMP
        PROPERTIES
            IMPORTED_LOCATION "${GMP_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}")
endif()

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY)
