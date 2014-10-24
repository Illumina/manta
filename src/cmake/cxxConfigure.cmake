#
# Manta
# Copyright (c) 2013-2014 Illumina, Inc.
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

include ("${THIS_MACROS_CMAKE}")

# Support for static linking
# Note that this implies that all libraries must be found with the
# exact file name (libXXX.a or libXXX.so)
#if    (THIS_FORCE_STATIC_LINK)
#    message(STATUS "All libraries will be statically linked")
#    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "-static")
#    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-static")
    # ensure that even if cmake decides to allow for dynamic libs resolution,
    # this gets overriden into static...
#    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS ${CMAKE_EXE_LINK_STATIC_CXX_FLAGS})
#    set(THIS_LIBRARY_PREFIX ${CMAKE_STATIC_LIBRARY_PREFIX})
#    set(THIS_LIBRARY_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
    # set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
#else  ()
#    set(THIS_LIBRARY_PREFIX "")
#    set(THIS_LIBRARY_SUFFIX "")
#endif ()

# optional support for gzip compression
static_find_library(ZLIB zlib.h z)
if    (HAVE_ZLIB)
    set  (THIS_ADDITIONAL_LIB ${THIS_ADDITIONAL_LIB} z)
    message(STATUS "gzip compression supported")
else  ()
    message(FATAL_ERROR "No support for gzip compression")
endif ()

# samtools 0.2.x forces pthreads in link:
find_package( Threads )


# Force static linking
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")

function(get_compiler_name_version compiler_name compiler_version)
    execute_process(COMMAND ${compiler_name} -dumpversion OUTPUT_VARIABLE this_version)
    STRING(REGEX REPLACE "(\r?\n)+$" "" this_version "${this_version}")
    set(${compiler_version} ${this_version} PARENT_SCOPE)
endfunction()

macro(get_compiler_version compiler_version)
    get_compiler_name_version(${CMAKE_CXX_COMPILER} compiler_version)
endmacro()

# clang doesn't make finding the version easy for us...
macro(get_clang_version compiler_version)
#    execute_process(COMMAND bash -c "${CMAKE_CXX_COMPILER} -v 2>&1 | awk '{printf $3; exit}'" OUTPUT_VARIABLE ${compiler_version})
    execute_process(COMMAND bash -c "echo | ${CMAKE_CXX_COMPILER} -dM -E - | awk '/__clang_version__/ {printf $3; exit}' | tr -d '\"'" OUTPUT_VARIABLE ${compiler_version})
endmacro()

macro(test_min_compiler compiler_version min_compiler_version compiler_label)
    if (${compiler_version} VERSION_LESS ${min_compiler_version})
        message (FATAL_ERROR "Unsupported version for ${compiler_label}: ${compiler_version}: "
                             "only versions >= ${min_compiler_version} are supported")
    endif ()
endmacro()


set(min_gxx_version "4.7")
set(min_clang_version "3.2")
set(min_intel_version "12.0") # guestimate based on intel support documentation

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    get_compiler_version(compiler_version)
    test_min_compiler(${compiler_version} "${min_gxx_version}" "g++")
    message (STATUS "using compiler: g++ version ${compiler_version}")

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    get_clang_version(compiler_version)
    test_min_compiler(${compiler_version} "${min_clang_version}" "clang++")
    message (STATUS "using compiler: clang++ version ${compiler_version}")

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    get_compiler_version(compiler_version)
    test_min_compiler(${compiler_version} "${min_intel_version}" "icpc") 
    message (STATUS "using compiler: Intel version ${compiler_version}")

    # for intel we also need to test the minimum version of g++ currently
    # in the path (because this is the stdc++ library that # intel will use):
    get_compiler_name_version("g++" gxx_compiler_version)
    test_min_compiler(${gxx_compiler_version} "${min_gxx_version}" "g++ libstdc++ (library used by icpc)")
    message (STATUS "using libstdc++: gnu version ${gxx_compiler_version}")

else ()
    message (STATUS "using compiler: ${CMAKE_CXX_COMPILER_ID}")
endif ()


#
# set compile flags, and modify by compiler/version:
#
set (GNU_COMPAT_COMPILER ( (CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")))

# start with warning flags:
if (GNU_COMPAT_COMPILER)
    set (CXX_WARN_FLAGS "-Wall -Wextra -Wshadow -Wunused -Wpointer-arith -Winit-self -pedantic -Wunused-parameter -Wundef -Wdisabled-optimization")
    if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        set (CXX_WARN_FLAGS "-Wredundant-decls")
    endif ()

    if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
        set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Wuninitialized")
    endif ()
endif ()

#
# add extra compiler specific flags:
#
if     (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if (NOT (${compiler_version} VERSION_LESS "4.2"))
        set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Wlogical-op")
    endif ()

    if (NOT (${compiler_version} VERSION_LESS "4.5"))
        # Force static linking of standard libraries:
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
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
        set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Wheader-guard -Wlogical-not-parentheses -Wloop-analysis")
    endif ()

    if (NOT (${compiler_version} VERSION_LESS "3.5"))
        set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Wextra-semi -Wdeprecated -Wmissing-variable-declarations")
    endif ()

    #### documentation of other possible warning flags from clang
    # set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Weverything -Wno-sign-conversion -Wno-weak-vtables -Wno-conversion -Wno-cast-align -Wno-padded -Wno-switch-enum -Wno-missing-noreturn -Wno-covered-switch-default -Wno-unreachable-code -Wno-global-constructors -Wno-exit-time-destructors")
    ### new disabled everything-warnings in clang 3.5 (or these might be c++11 warnings):
    # set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Wno-c++98-compat -Wno-documentation-unknown-command -Wno-old-style-cast -Wno-unused-member-function -Wno-documentation -Wno-float-equal")

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # suppress errors in boost headers:
    set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -diag-disable 177,193,869,1599,3280")

    set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Wunused-variable -Wpointer-arith -Wuninitialized")

    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
    #set (CXX_WARN_FLAGS "${CXX_WARN_FLAGS} -Wmissing-prototypes -Wmissing-declarations -Wunused-variable -Wpointer-arith -Wuninitialized")
endif()


# The NDEBUG macro is intentionally removed from release. One discussion on this is:
# http://www.drdobbs.com/an-exception-or-a-bug/184401686

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_WARN_FLAGS}")

if (GNU_COMPAT_COMPILER)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
    set (CMAKE_CXX_FLAGS_RELEASE "-O3 -fomit-frame-pointer")
    set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
    set (CMAKE_CXX_FLAGS_ASAN "-O1 -g -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls")
    #set (CMAKE_CXX_FLAGS_PROFILE "-O0 -g -pg -fprofile-arcs -ftest-coverage")
endif()

# if ASan build type is requested, check that the compiler supports it:
if (CMAKE_BUILD_TYPE STREQUAL "ASan")
    set (IS_ASAN_SUPPORTED false)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        if (NOT (${compiler_version} VERSION_LESS "4.8"))
            set (IS_ASAN_SUPPORTED true)
        endif ()
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        if (NOT (${compiler_version} VERSION_LESS "3.1"))
            set (IS_ASAN_SUPPORTED true)
        endif ()
    endif ()

    if (NOT ${IS_ASAN_SUPPORTED})
        message(FATAL_ERROR "Address sanitizer build type requested, but this is not supported by compiler.")
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
        message (STATUS "building in developer mode: treating compiler warnings as errors")
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
    ## Prevent using 80bit registers (more consistent rounding)
    ##
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffloat-store")
  endif ()

endif()

set(THIS_CXX_CONFIG_H_DIR ${CMAKE_CURRENT_BINARY_DIR}/lib)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lib/common/config.h.in ${THIS_CXX_CONFIG_H_DIR}/common/config.h @ONLY)


#
# include dirs:
#
set (THIS_CXX_BEFORE_SYSTEM_INCLUDES "${Boost_INCLUDE_DIRS}" "${HTSLIB_DIR}" "${SAMTOOLS_DIR}")
set (THIS_CXX_ALL_INCLUDES "${CMAKE_SOURCE_DIR}/c++/lib")
