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
## Top level cmake file
##
## author Come Raczy
##
################################################################################

cmake_minimum_required(VERSION 2.8.0)

message (STATUS "**** Initializing project cmake configuration ****")

enable_testing()

#
# setup project/version
#

#project (MANTA)
set(MANTA_PROJECT_NAME "manta")

# get version
set(SOURCE_BOOTSTRAP_DIR "${CMAKE_SOURCE_DIR}/cmake/bootstrap")

execute_process(COMMAND "${SOURCE_BOOTSTRAP_DIR}/getVersion.bash" OUTPUT_VARIABLE MANTA_VERSION)
set(MANTA_VERSION_FULL "${MANTA_PROJECT_NAME}-${MANTA_VERSION}")
message (STATUS "MANTA_VERSION: ${MANTA_VERSION}")


#
# setup build type
#
if (NOT CMAKE_BUILD_TYPE)
    set(DEFAULT_CMAKE_BUILD_TYPE RelWithDebInfo)
    set(CMAKE_BUILD_TYPE ${DEFAULT_CMAKE_BUILD_TYPE} CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo (default: ${DEFAULT_CMAKE_BUILD_TYPE})" FORCE)
endif()
message (STATUS "BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

#
# setup cpack
#
if (CPACK_GENERATOR)
    message (STATUS "Configuring to produce the following package types: ${CPACK_GENERATOR}")
    SET(CPACK_PACKAGE_VENDOR "Illumina")
    SET(CPACK_PACKAGE_VERSION "${MANTA_VERSION}")
    SET(CPACK_SET_DESTDIR ON)
    INCLUDE(CPack)
endif ()

# Installation directories
if    (NOT CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX "/usr/local")
endif ()

# check that the installation prefix is acceptable
function(check_prefix)
    string(SUBSTRING "${CMAKE_INSTALL_PREFIX}" 0 1 slash)
    string(COMPARE NOTEQUAL "/" "${slash}" isslash)
    if    (isslash)
        install(CODE "
            string(LENGTH \"\$ENV{DESTDIR}\" ddlen)
            if (ddlen)
                message (FATAL_ERROR \"Installation prefix must begin with '/' if DESTDIR is set.\")
            endif ()
        ")
    endif ()
endfunction()

check_prefix()

set(MANTA_EXEC_PREFIX ".")
set(MANTA_BINDIR "${MANTA_EXEC_PREFIX}/bin")
set(MANTA_LIBDIR "${MANTA_EXEC_PREFIX}/lib/${MANTA_VERSION_FULL}")
set(MANTA_PYTHON_LIBDIR "${MANTA_LIBDIR}/python")
set(MANTA_LIBEXECDIR "${MANTA_EXEC_PREFIX}/libexec/${MANTA_VERSION_FULL}")
set(MANTA_DATADIR "share/${MANTA_VERSION_FULL}")
set(MANTA_ETCDIR "${MANTA_EXEC_PREFIX}/etc/${MANTA_VERSION_FULL}")
set(MANTA_DOCDIR "doc/${MANTA_VERSION_FULL}")

message (STATUS "install prefix: ${CMAKE_INSTALL_PREFIX}")

#add_custom_target(MANTA_CXX)
add_custom_target(MANTA_OPT ALL)

set(BOOST_REDIST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../redist")

# required boost libraries
set (MANTA_BOOST_VERSION 1.49.0)
set (MANTA_BOOST_COMPONENTS date_time filesystem iostreams program_options
                            regex system unit_test_framework)

# the name given to boost.build and the library name are the same for all libraries, except
# for test, so we need two lists now:
set (MANTA_BOOST_BUILD_COMPONENTS date_time filesystem iostreams program_options
                                 regex system test)
set (Boost_USE_MULTITHREADED OFF)
set (Boost_USE_STATIC_LIBS ON)

set (MANTA_GLOBALS_CMAKE "${CMAKE_SOURCE_DIR}/cmake/globals.cmake")
set (MANTA_MACROS_CMAKE "${CMAKE_SOURCE_DIR}/cmake/macros.cmake")

# redist includes
include ("${CMAKE_SOURCE_DIR}/cmake/boost.cmake")

# force unit tests at compile time, and set them to be more verbose than
# cmake's default "test" target:
#add_custom_target(MANTA_UNITTESTS
#    ALL
#    COMMAND ${CMAKE_CTEST_COMMAND} --force-new-ctest-process --output-on-failure
#    )

add_subdirectory (cmake/preInstall)
add_subdirectory (opt)
add_subdirectory (c++)
#add_subdirectory (python)
add_subdirectory (cmake/postInstall)



