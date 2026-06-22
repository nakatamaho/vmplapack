# Copyright (c) 2026, NAKATA Maho
# SPDX-License-Identifier: BSD-2-Clause

cmake_minimum_required(VERSION 3.16)

project(vmplapack
    VERSION 0.1.0
    LANGUAGES CXX)

option(VMPLAPACK_ENABLE_MPFR "Enable the MPFR tier." ON)
option(VMPLAPACK_ENABLE_TESTS "Build vMPLAPACK tests." ON)
option(VMPLAPACK_ENABLE_EXAMPLES "Build vMPLAPACK examples." ON)
option(VMPLAPACK_USE_DEKKER_TWOPRODUCT "Use Dekker TwoProduct instead of FMA." OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(RoundingControl)
include(CheckCXXSourceCompiles)

add_library(vmplapack INTERFACE)
add_library(vmplapack::vmplapack ALIAS vmplapack)

target_compile_features(vmplapack INTERFACE cxx_std_17)
target_include_directories(vmplapack
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
target_link_libraries(vmplapack INTERFACE vmplapack::rounding_control)

if(VMPLAPACK_ENABLE_MPFR)
    set(GMPFRXX_MKII_INCLUDE_DIR "" CACHE PATH "Path to the gmpfrxx_mkII include directory")
    if(NOT GMPFRXX_MKII_INCLUDE_DIR)
        find_path(VMPLAPACK_DETECTED_GMPFRXX_MKII_INCLUDE_DIR
            NAMES mpfrxx_mkII.h gmpfrxx_mkII.h
            PATHS
                "${CMAKE_CURRENT_SOURCE_DIR}/extern/gmpfrxx_mkII/include"
                "/home/docker/gmpfrxx_mkII/include")
        if(VMPLAPACK_DETECTED_GMPFRXX_MKII_INCLUDE_DIR)
            set(GMPFRXX_MKII_INCLUDE_DIR "${VMPLAPACK_DETECTED_GMPFRXX_MKII_INCLUDE_DIR}"
                CACHE PATH "Path to the gmpfrxx_mkII include directory" FORCE)
        endif()
    endif()
    if(NOT GMPFRXX_MKII_INCLUDE_DIR)
        message(FATAL_ERROR "M2 MPFR API probe failed: missing gmpfrxx_mkII headers (set GMPFRXX_MKII_INCLUDE_DIR)")
    endif()

    find_path(GMP_INCLUDE_DIR
        NAMES gmp.h
        PATHS
            "/home/docker/gmpfrxx_mkII/build_release/_deps/gmpfrxx_mkii/GMP/include"
            "/home/docker/mplapack/external/i/GMP/include")
    find_path(MPFR_INCLUDE_DIR
        NAMES mpfr.h
        PATHS
            "/home/docker/gmpfrxx_mkII/build_release/_deps/gmpfrxx_mkii/MPFR/include"
            "/home/docker/mplapack/external/i/MPFR/include")
    find_library(GMP_LIBRARY
        NAMES gmp
        PATHS
            "/home/docker/gmpfrxx_mkII/build_release/_deps/gmpfrxx_mkii/GMP/lib"
            "/home/docker/mplapack/external/i/GMP/lib")
    find_library(MPFR_LIBRARY
        NAMES mpfr
        PATHS
            "/home/docker/gmpfrxx_mkII/build_release/_deps/gmpfrxx_mkii/MPFR/lib"
            "/home/docker/mplapack/external/i/MPFR/lib")
    find_path(MPC_INCLUDE_DIR
        NAMES mpc.h
        PATHS
            "/home/docker/gmpfrxx_mkII/build_release/_deps/gmpfrxx_mkii/MPC/include"
            "/home/docker/mplapack/external/i/MPC/include")
    find_library(MPC_LIBRARY
        NAMES mpc
        PATHS
            "/home/docker/gmpfrxx_mkII/build_release/_deps/gmpfrxx_mkii/MPC/lib"
            "/home/docker/mplapack/external/i/MPC/lib")
    if(NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARY)
        message(FATAL_ERROR "M2 MPFR API probe failed: missing GMP headers or library")
    endif()
    if(NOT MPFR_INCLUDE_DIR OR NOT MPFR_LIBRARY)
        message(FATAL_ERROR "M2 MPFR API probe failed: missing MPFR headers or library")
    endif()

    add_library(vmplapack_mpfr_deps INTERFACE)
    add_library(vmplapack::mpfr_deps ALIAS vmplapack_mpfr_deps)
    target_include_directories(vmplapack_mpfr_deps
        INTERFACE
            "${GMPFRXX_MKII_INCLUDE_DIR}"
            "${GMP_INCLUDE_DIR}"
            "${MPFR_INCLUDE_DIR}")
    target_link_libraries(vmplapack_mpfr_deps INTERFACE "${MPFR_LIBRARY}" "${GMP_LIBRARY}")

    target_compile_definitions(vmplapack INTERFACE VMPLAPACK_ENABLE_MPFR)
    target_link_libraries(vmplapack INTERFACE vmplapack::mpfr_deps)

    set(_vmplapack_probe_required_includes "${GMPFRXX_MKII_INCLUDE_DIR};${GMP_INCLUDE_DIR};${MPFR_INCLUDE_DIR}")
    set(_vmplapack_probe_required_libraries "${MPFR_LIBRARY};${GMP_LIBRARY}")
    if(MSVC)
        set(_vmplapack_probe_required_flags "")
    else()
        set(_vmplapack_probe_required_flags "-std=c++17")
    endif()

    function(vmplapack_probe_mpfr_api symbol source)
        string(MAKE_C_IDENTIFIER "${symbol}" probe_id)
        set(CMAKE_REQUIRED_INCLUDES "${_vmplapack_probe_required_includes}")
        set(CMAKE_REQUIRED_LIBRARIES "${_vmplapack_probe_required_libraries}")
        set(CMAKE_REQUIRED_FLAGS "${_vmplapack_probe_required_flags}")
        set(CMAKE_REQUIRED_QUIET TRUE)
        check_cxx_source_compiles("${source}" "VMPLAPACK_MPFR_PROBE_${probe_id}")
        if(NOT "${VMPLAPACK_MPFR_PROBE_${probe_id}}")
            message(FATAL_ERROR "M2 MPFR API probe failed: missing symbol ${symbol}")
        endif()
    endfunction()

    vmplapack_probe_mpfr_api("mpfrxx::mpfr_class" [[
        #include <mpfrxx_mkII.h>
        int main() {
            mpfrxx::mpfr_class x;
            (void)x;
            return 0;
        }
    ]])
    vmplapack_probe_mpfr_api("mpfrxx::fma" [[
        #include <mpfrxx_mkII.h>
        int main() {
            mpfrxx::mpfr_class a(1.0);
            mpfrxx::mpfr_class b(2.0);
            mpfrxx::mpfr_class c(3.0);
            mpfrxx::mpfr_class y = mpfrxx::fma(a, b, c);
            return mpfr_number_p(y.mpfr_data()) ? 0 : 1;
        }
    ]])
    vmplapack_probe_mpfr_api("mpfrxx::rounding_mode_scope" [[
        #include <mpfrxx_mkII.h>
        int main() {
            mpfrxx::rounding_mode_scope up(MPFR_RNDU);
            mpfrxx::rounding_mode_scope down(MPFR_RNDD);
            return 0;
        }
    ]])
    vmplapack_probe_mpfr_api("mpfr_class::mpfr_data/get_mpfr_t" [[
        #include <mpfrxx_mkII.h>
        int main() {
            mpfrxx::mpfr_class x(1.0);
            mpfr_srcptr a = x.mpfr_data();
            mpfr_srcptr b = x.get_mpfr_t();
            return (a != nullptr && b != nullptr) ? 0 : 1;
        }
    ]])
    vmplapack_probe_mpfr_api("mpfr_class::with_precision" [[
        #include <mpfrxx_mkII.h>
        int main() {
            mpfrxx::mpfr_class x = mpfrxx::mpfr_class::with_precision(128);
            mpfrxx::mpfr_class y = mpfrxx::mpfr_class::with_precision(128, 1.0);
            return (mpfr_get_prec(x.mpfr_data()) == 128 && mpfr_get_prec(y.mpfr_data()) == 128) ? 0 : 1;
        }
    ]])
    vmplapack_probe_mpfr_api("mpfrxx::default_prec/set_default_precision_bits" [[
        #include <mpfrxx_mkII.h>
        int main() {
            mpfrxx::set_default_precision_bits(192);
            return mpfrxx::default_prec() == 192 ? 0 : 1;
        }
    ]])
endif()

if(VMPLAPACK_USE_DEKKER_TWOPRODUCT)
    target_compile_definitions(vmplapack INTERFACE VMPLAPACK_USE_DEKKER_TWOPRODUCT)
endif()

if(VMPLAPACK_ENABLE_TESTS)
    include(CTest)
    if(BUILD_TESTING)
        add_executable(test_Rarith_environment tests/test_Rarith_environment.cpp)
        target_link_libraries(test_Rarith_environment PRIVATE vmplapack::vmplapack)
        add_test(NAME test_Rarith_environment COMMAND test_Rarith_environment)

        add_executable(test_Rarith_eft tests/test_Rarith_eft.cpp)
        target_link_libraries(test_Rarith_eft PRIVATE vmplapack::vmplapack)
        add_test(NAME test_Rarith_eft COMMAND test_Rarith_eft)
        if(VMPLAPACK_ENABLE_MPFR)
            add_test(NAME test_Rarith_eft_mpfr_w53 COMMAND test_Rarith_eft)
            set_tests_properties(test_Rarith_eft_mpfr_w53
                PROPERTIES ENVIRONMENT "MPFRXX_DEFAULT_PRECISION_BITS=53")
            add_test(NAME test_Rarith_eft_mpfr_w512 COMMAND test_Rarith_eft)
            set_tests_properties(test_Rarith_eft_mpfr_w512
                PROPERTIES ENVIRONMENT "MPFRXX_DEFAULT_PRECISION_BITS=512")

            add_executable(test_Rdot_oracle_generators tests/test_Rdot_oracle_generators.cpp)
            target_link_libraries(test_Rdot_oracle_generators PRIVATE vmplapack::vmplapack)
            add_test(NAME test_Rdot_oracle_generators COMMAND test_Rdot_oracle_generators)
        endif()
    endif()
endif()

if(VMPLAPACK_ENABLE_EXAMPLES)
    add_executable(example_m0_midrad examples/example_m0_midrad.cpp)
    target_link_libraries(example_m0_midrad PRIVATE vmplapack::vmplapack)

    add_executable(example_m1_twoproduct_residual examples/example_m1_twoproduct_residual.cpp)
    target_link_libraries(example_m1_twoproduct_residual PRIVATE vmplapack::vmplapack)

    add_executable(example_m1_high_condition_dot_preview examples/example_m1_high_condition_dot_preview.cpp)
    target_link_libraries(example_m1_high_condition_dot_preview PRIVATE vmplapack::vmplapack)
endif()
