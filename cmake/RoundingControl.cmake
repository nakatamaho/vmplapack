# Copyright (c) 2026, NAKATA Maho
# SPDX-License-Identifier: BSD-2-Clause

add_library(vmplapack_rounding_control INTERFACE)
add_library(vmplapack::rounding_control ALIAS vmplapack_rounding_control)

if(MSVC)
    target_compile_options(vmplapack_rounding_control
        INTERFACE
            /W4
            /WX
            /fp:strict)
else()
    target_compile_options(vmplapack_rounding_control
        INTERFACE
            -Wall
            -Wextra
            -Werror
            -fno-fast-math
            -frounding-math
            -ffp-contract=off)

    if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(i[3-6]86|x86)$")
        target_compile_options(vmplapack_rounding_control
            INTERFACE
                -mfpmath=sse
                -msse2)
    endif()
endif()
