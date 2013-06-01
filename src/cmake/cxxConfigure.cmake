#
# Manta
# Copyright (c) 2013 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/downloads/sequencing/licenses/>.
#

################################################################################
##
## CMake configuration file for c++ executables
##
## author Come Raczy
##
################################################################################


INCLUDE(TestBigEndian)
TEST_BIG_ENDIAN(IS_BIG_ENDIAN)

INCLUDE(CheckFunctionExists)

find_path(HAVE_INTTYPES_H  inttypes.h)
find_path(HAVE_MEMORY_H    memory.h)
find_path(HAVE_STDINT_H    stdint.h)
find_path(HAVE_STDLIB_H    stdlib.h)
find_path(HAVE_STRING_H    string.h)
find_path(HAVE_STRINGS_H   strings.h)
find_path(HAVE_UNISTD_H    unistd.h)
find_path(HAVE_SYS_STAT_H  sys/stat.h)
find_path(HAVE_SYS_TYPES_H sys/types.h)

set (CMAKE_REQUIRED_LIBRARIES m)
check_function_exists(floorf HAVE_FLOORF)
check_function_exists(round  HAVE_ROUND)
check_function_exists(roundf HAVE_ROUNDF)
check_function_exists(powf HAVE_POWF)
check_function_exists(erf HAVE_ERF)
check_function_exists(erf HAVE_ERFF)
check_function_exists(erfc HAVE_ERFC)
check_function_exists(erfc HAVE_ERFCF)

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

# optional support for bzip2 compression
static_find_library(BZIP2 bzlib.h bz2)
if    (HAVE_BZIP2)
    set(HAVE_BZLIB HAVE_BZIP2)
    set(MANTA_ADDITIONAL_LIB ${MANTA_ADDITIONAL_LIB} bz2)
    message(STATUS "bzip2 compression supported")
else  ()
    message(FATAL_ERROR "No support for bzip2 compression")
endif ()

static_find_boost(${MANTA_BOOST_VERSION} "${MANTA_BOOST_COMPONENTS}")

static_find_library(CPGPLOT cpgplot.h cpgplot)
static_find_library(PGPLOT cpgplot.h pgplot)
static_find_library(X11 X.h X11)

set (CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} -Wall -Wextra -Wunused -Wno-long-long -Wsign-compare -Wpointer-arith" CACHE STRING "g++ flags" FORCE)
set (CMAKE_CXX_FLAGS_DEBUG "-O0 -DNDEBUG -g " CACHE STRING "g++ flags" FORCE)
set (CMAKE_CXX_FLAGS_PROFILE "-O0 -DNDEBUG -g -pg -fprofile-arcs -ftest-coverage" CACHE STRING "g++ flags" FORCE)
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "g++ flags" FORCE)
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g" CACHE STRING "g++ flags" )

# Force static linking
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE version)
    
    string(REGEX REPLACE "^([0-9])\\.[0-9]\\.[0-9]" "\\1" major_version ${version})
    string(REGEX REPLACE "^[0-9]\\.([0-9])\\.[0-9]" "\\1" minor_version ${version})
    string(REGEX REPLACE "^[0-9]\\.[0-9]\\.([0-9])" "\\1" patch_version ${version})
    if    (major_version LESS 4 OR major_version EQUAL 4 AND minor_version LESS 1)
        message (FATAL_ERROR "Unsupported GNU C++ compiler: g++ version ${version}: "
                             "only g++ versions >= 4.1.0 are supported")
    endif ()

    set("${CMAKE_CXX_COMPILER_ID}${major_version}" true)
    set("${CMAKE_CXX_COMPILER_ID}${major_version}${minor_version}" true)
    set("${CMAKE_CXX_COMPILER_ID}${major_version}${minor_version}${patch_version}" true)
    message (STATUS "using compiler: gcc version ${version}")
endif ()

##
## Suppress spurious warnings in less recent compilers
##
if    (NOT GNU42)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter ")
endif ()

if    (GNU412 OR GNU42 OR GNU43)
    ## Before 4.1.2, pedantic breaks on boost lambda expressions
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic ")
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

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lib/common/config.h.in ${MANTA_CXX_CONFIG_H_DIR}/config.h @ONLY)
