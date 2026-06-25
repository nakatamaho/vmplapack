# Copyright (c) 2026, NAKATA Maho
# SPDX-License-Identifier: BSD-2-Clause

if(NOT DEFINED VMPLAPACK_CONSUMER_BUILD_DIR)
    message(FATAL_ERROR "VMPLAPACK_CONSUMER_BUILD_DIR is required")
endif()

set(_compile_commands "${VMPLAPACK_CONSUMER_BUILD_DIR}/compile_commands.json")
if(NOT EXISTS "${_compile_commands}")
    message(FATAL_ERROR "compile_commands.json was not generated")
endif()

file(READ "${_compile_commands}" _commands)
string(FIND "${_commands}" "consumer_verified_dot.cpp" _source_pos)
if(_source_pos EQUAL -1)
    message(FATAL_ERROR "consumer_verified_dot.cpp was not found in compile_commands.json")
endif()

if(VMPLAPACK_CONSUMER_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    foreach(_flag IN ITEMS -fno-fast-math -frounding-math -ffp-contract=off)
        string(FIND "${_commands}" "${_flag}" _flag_pos)
        if(_flag_pos EQUAL -1)
            message(FATAL_ERROR "required inherited FP flag was not found: ${_flag}")
        endif()
    endforeach()
elseif(VMPLAPACK_CONSUMER_COMPILER_ID STREQUAL "MSVC")
    string(FIND "${_commands}" "/fp:strict" _fp_pos)
    if(_fp_pos EQUAL -1)
        message(FATAL_ERROR "required inherited MSVC FP flag was not found: /fp:strict")
    endif()
else()
    message(FATAL_ERROR "unsupported consumer compiler ID: ${VMPLAPACK_CONSUMER_COMPILER_ID}")
endif()
