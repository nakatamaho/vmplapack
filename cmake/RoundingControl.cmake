# Copyright (c) 2026, NAKATA Maho
# SPDX-License-Identifier: BSD-2-Clause

add_library(vmplapack_rounding_control INTERFACE)
add_library(vmplapack::rounding_control ALIAS vmplapack_rounding_control)
set_target_properties(vmplapack_rounding_control PROPERTIES EXPORT_NAME rounding_control)

if(NOT MSVC AND NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    message(FATAL_ERROR "vMPLAPACK requires strict floating-point flags for this compiler; unsupported compiler ID: ${CMAKE_CXX_COMPILER_ID}")
endif()

set(_vmplapack_x86_32 FALSE)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(i[3-6]86|x86)$")
    set(_vmplapack_x86_32 TRUE)
endif()

target_compile_options(vmplapack_rounding_control
    INTERFACE
        "$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/W4;/WX;/fp:strict>"
        "$<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>:-Wall;-Wextra;-Werror;-fno-fast-math;-frounding-math;-ffp-contract=off>"
        "$<$<AND:$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>,$<BOOL:${_vmplapack_x86_32}>>:-mfpmath=sse;-msse2>")
