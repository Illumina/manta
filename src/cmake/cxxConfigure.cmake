#
# Manta
# Copyright (c) 2013 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
#

################################################################################
##
## CMake configuration file for c++ executables
##
## author Come Raczy
##
################################################################################


INCLUDE(CheckFunctionExists)

find_path(HAVE_INTTYPES_H  inttypes.h)
find_path(HAVE_MEMORY_H    memory.h)
find_path(HAVE_STDINT_H    stdint.h)
find_path(HAVE_STDLIB_H    stdlib.h)
find_path(HAVE_STRING_H    string.h)
find_path(HAVE_STRINGS_H   strings.h)
find_path(HAVE_UNISTD_H    unistd.h)

set (CMAKE_REQUIRED_LIBRARIES m)
check_function_exists(floorf HAVE_FLOORF)
check_function_exists(round  HAVE_ROUND)
check_function_exists(roundf HAVE_ROUNDF)
check_function_exists(powf HAVE_POWF)

include ("${MANTA_MACROS_CMAKE}")

# Support for static linking
# Note that this implies that all libraries must be found with the
# exact file name (libXXX.a or libXXX.so)
if    (MANTA_FORCE_STATIC_LINK)
    message(STATUS "All libraries will be statically linked")
    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "-static")
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-static")
    # ensure that even if cmake decides to allow for dynamic libs resolution,
    # this gets overriden into static...
    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS ${CMAKE_EXE_LINK_STATIC_CXX_FLAGS})
    set(MANTA_LIBRARY_PREFIX ${CMAKE_STATIC_LIBRARY_PREFIX})
    set(MANTA_LIBRARY_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
else  ()
    set(MANTA_LIBRARY_PREFIX "")
    set(MANTA_LIBRARY_SUFFIX "")
endif ()

# optional support for gzip compression
static_find_library(ZLIB zlib.h z)
if    (HAVE_ZLIB)
    set  (MANTA_ADDITIONAL_LIB ${MANTA_ADDITIONAL_LIB} z)
    message(STATUS "gzip compression supported")
else  ()
    message(FATAL_ERROR "No support for gzip compression")
endif ()

static_find_boost(${MANTA_BOOST_VERSION} "${MANTA_BOOST_COMPONENTS}")

# Force static linking
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")

macro(get_compiler_version compiler_version)
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE ${compiler_version})
    STRING(REGEX REPLACE "(\r?\n)+$" "" ${compiler_version} "${${compiler_version}}")
endmacro()

# clang doesn't make finding the version easy for us...
macro(get_clang_version compiler_version)
    execute_process(COMMAND bash -c "${CMAKE_CXX_COMPILER} -v 2>&1 | awk '{printf $3; exit}'" OUTPUT_VARIABLE ${compiler_version})
endmacro()

macro(test_min_compiler compiler_version min_compiler_version compiler_label)
    if (${compiler_version} VERSION_LESS ${min_compiler_version})
        message (FATAL_ERROR "Unsupported ${compiler_label} version: ${compiler_version}: "
                             "only versions >= ${min_compiler_version} are supported")
    endif ()
endmacro()


if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    get_compiler_version(compiler_version)
    test_min_compiler(${compiler_version} "4.1.2" "g++")
    message (STATUS "using compiler: g++ version ${compiler_version}")

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    get_clang_version(compiler_version)
    test_min_compiler(${compiler_version} "3.2" "clang++")
    message (STATUS "using compiler: clang++ version ${compiler_version}")
else ()
    message (STATUS "using compiler: ${CMAKE_CXX_COMPILER_ID}")
endif ()


#
# set compile flags, and modify by compiler/version:
#
set (GNU_COMPAT_COMPILER ( (CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang") ))

# start with warning flags:
if (GNU_COMPAT_COMPILER)
    set (CXX_WARN_FLAGS "-Wall -Wextra -Wshadow -Wunused -Wpointer-arith -Winit-self -Wredundant-decls -pedantic -Wunused-parameter -Wundef -Wdisabled-optimization")

    if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
        set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Wuninitialized")
    endif ()
endif ()

#
# add extra compiler specific warnings:
#
if     (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if (NOT (${compiler_version} VERSION_LESS "4.2"))
        set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Wlogical-op")
    endif ()

    if (NOT ((${compiler_version} VERSION_LESS "4.7") OR (${compiler_version} VERSION_GREATER "4.7")))
        # switching off warning about unused function because otherwise compilation will fail with g++ 4.7.3 in Ubuntu
        set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Wno-unused-function")
    endif ()

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Wmissing-prototypes -Wunused-exception-parameter -Wbool-conversion -Wempty-body -Wimplicit-fallthrough -Wsizeof-array-argument -Wstring-conversion")

    if (NOT (${compiler_version} VERSION_LESS "3.3"))
        set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Woverloaded-shift-op-parentheses")
    endif ()

    if (NOT (${compiler_version} VERSION_LESS "3.4"))
        # wait for 3.4 to be released before turning these on by default
        #set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Wheader-guard -Wlogical-not-parentheses -Wloop-analysis -Wunique-enum")
    endif ()

    # documentation of other possible warning flags from clang
    #set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Weverything -Wno-sign-conversion -Wno-weak-vtables -Wno-conversion -Wno-cast-align -Wno-padded -Wno-switch-enum -Wno-missing-noreturn -Wno-covered-switch-default -Wno-unreachable-code -Wno-global-constructors -Wno-exit-time-destructors")
endif()


# The NDEBUG macro is intentionally removed from release. One discussion on this is:
# http://www.drdobbs.com/an-exception-or-a-bug/184401686

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_WARN_FLAGS}")

if (GNU_COMPAT_COMPILER)
    set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
    set (CMAKE_CXX_FLAGS_RELEASE "-O3 -fomit-frame-pointer")
    set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
    #set (CMAKE_CXX_FLAGS_PROFILE "-O0 -g -pg -fprofile-arcs -ftest-coverage")
endif()


# add address sanitizer to debug mode:
set (USE_ADDRESS_SANITIZER false) # if true, turn on Address Sanitizer in debug for compilers which support this:

if (${USE_ADDRESS_SANITIZER})
    set (IS_ASAN false)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        if (NOT (${compiler_version} VERSION_LESS "4.8"))
            set (IS_ASAN true)
        endif ()
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        if (NOT (${compiler_version} VERSION_LESS "3.1"))
            set (IS_ASAN true)
        endif ()
    endif ()

    if (${IS_ASAN})
        set (CMAKE_CXX_FLAGS_DEBUG "-fsanitize=address -fno-omit-frame-pointer ${CMAKE_CXX_FLAGS_DEBUG}")
    endif ()
endif ()


if (GNU_COMPAT_COMPILER)

  if (${DEVELOPER_MODE})
    # some compiler versions will produce warnings with no reasonable workaround,
    # turn Werror off in this case
    #
    # a very common example are warnings from boost generated despite this library
    # being identified as a system header
    #
    set(IS_WERROR true)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        if (${compiler_version} VERSION_LESS "4.2")
            set(IS_WERROR false)
        endif ()
    endif ()

    if(${IS_WERROR})
        message (STATUS "building in developer mode: treating compiler errors as warnings")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
    endif ()
  endif ()

  if (CMAKE_SYSTEM_PROCESSOR MATCHES "^i[67]86$")
    ##
    ## Use scalar floating point instructions from the SSE instruction set.
    ## Note: Pentium3 SSE supports only single precision arithmetics
    ##
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse -mfpmath=sse")
  elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "^i[345]86$")
    ##
    ## Prevent using 80bits registers (more consistent rounding)
    ##
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffloat-store")
  endif ()

endif()

set(MANTA_CXX_CONFIG_H_DIR ${CMAKE_CURRENT_BINARY_DIR}/lib/common)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lib/common/config.h.in ${MANTA_CXX_CONFIG_H_DIR}/config.h @ONLY)


#
# include dirs:
#
set (MANTA_CXX_BEFORE_SYSTEM_INCLUDES "${Boost_INCLUDE_DIRS}" "${SAMTOOLS_DIR}" "${TABIX_DIR}")
set (MANTA_CXX_ALL_INCLUDES "${CMAKE_SOURCE_DIR}/c++/lib")
