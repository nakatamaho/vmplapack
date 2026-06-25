# Copyright (c) 2026, NAKATA Maho
# SPDX-License-Identifier: BSD-2-Clause

cmake_minimum_required(VERSION 3.16)

project(vmplapack_consumer
    LANGUAGES CXX)

option(VMPLAPACK_CONSUMER_EXPECT_MPFR "Require the installed package to expose the MPFR tier." OFF)
option(VMPLAPACK_CONSUMER_FORCE_FAST_MATH "Add -ffast-math to verify the umbrella guard rejects it." OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(vmplapack 0.1.0 EXACT REQUIRED CONFIG)

function(vmplapack_consumer_append_plain_target_property output_variable target_name property_name flag_prefix)
    if(NOT TARGET "${target_name}")
        return()
    endif()

    get_target_property(_vmplapack_consumer_values "${target_name}" "${property_name}")
    if(NOT _vmplapack_consumer_values)
        return()
    endif()

    set(_vmplapack_consumer_result "${${output_variable}}")
    foreach(_vmplapack_consumer_value IN LISTS _vmplapack_consumer_values)
        if(_vmplapack_consumer_value MATCHES "^\\$<")
            continue()
        endif()
        list(APPEND _vmplapack_consumer_result "${flag_prefix}${_vmplapack_consumer_value}")
    endforeach()
    set("${output_variable}" "${_vmplapack_consumer_result}" PARENT_SCOPE)
endfunction()

if(VMPLAPACK_CONSUMER_FORCE_FAST_MATH)
    if(MSVC)
        message(FATAL_ERROR "vMPLAPACK correctly rejects /fp:fast consumers through the strict FP target policy")
    endif()

    set(_vmplapack_fast_math_compile_flags -ffast-math)
    foreach(_vmplapack_target IN ITEMS
            vmplapack::vmplapack
            vmplapack::mpfr_deps
            MPFR::MPFR
            GMP::GMP)
        vmplapack_consumer_append_plain_target_property(
            _vmplapack_fast_math_compile_flags
            "${_vmplapack_target}"
            INTERFACE_INCLUDE_DIRECTORIES
            "-I")
        vmplapack_consumer_append_plain_target_property(
            _vmplapack_fast_math_compile_flags
            "${_vmplapack_target}"
            INTERFACE_COMPILE_DEFINITIONS
            "-D")
    endforeach()

    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fast_math_negative")
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/fast_math_negative/fast_math_consumer.cpp" [[
#include <vmplapack/vmplapack.h>
int main() { return 0; }
]])
    try_compile(VMPLAPACK_CONSUMER_FAST_MATH_COMPILED
        "${CMAKE_CURRENT_BINARY_DIR}/fast_math_negative/build"
        "${CMAKE_CURRENT_BINARY_DIR}/fast_math_negative/fast_math_consumer.cpp"
        CMAKE_FLAGS
            "-DCMAKE_CXX_STANDARD:STRING=17"
            "-DCMAKE_CXX_STANDARD_REQUIRED:BOOL=ON"
        COMPILE_DEFINITIONS
            ${_vmplapack_fast_math_compile_flags}
        OUTPUT_VARIABLE VMPLAPACK_CONSUMER_FAST_MATH_OUTPUT)
    if(VMPLAPACK_CONSUMER_FAST_MATH_COMPILED)
        message(FATAL_ERROR "consumer -ffast-math compile was not rejected by include/vmplapack/vmplapack.h")
    endif()
    if(NOT VMPLAPACK_CONSUMER_FAST_MATH_OUTPUT MATCHES "vMPLAPACK requires IEEE-754 semantics")
        message(FATAL_ERROR "consumer -ffast-math compile failed, but not through the vMPLAPACK __FAST_MATH__ guard")
    endif()
    message(FATAL_ERROR "vMPLAPACK correctly rejected a consumer that added -ffast-math")
endif()

add_executable(consumer_verified_dot consumer_verified_dot.cpp)
target_link_libraries(consumer_verified_dot PRIVATE vmplapack::vmplapack)

if(VMPLAPACK_CONSUMER_EXPECT_MPFR)
    target_compile_definitions(consumer_verified_dot PRIVATE VMPLAPACK_CONSUMER_EXPECT_MPFR)
endif()


enable_testing()
add_test(NAME consumer_verified_dot COMMAND consumer_verified_dot)
add_test(NAME consumer_fp_flags
    COMMAND "${CMAKE_COMMAND}"
        -DVMPLAPACK_CONSUMER_BUILD_DIR=${CMAKE_BINARY_DIR}
        -DVMPLAPACK_CONSUMER_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}
        -P "${CMAKE_CURRENT_SOURCE_DIR}/verify_fp_flags.cmake")
