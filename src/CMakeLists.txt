#
# Manta
# Copyright (c) 2013-2015 Illumina, Inc.
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
## Top level cmake file
##
## author Come Raczy
##
################################################################################

cmake_minimum_required(VERSION 2.8.0)

message (STATUS "==== Initializing project cmake configuration ====")

enable_testing()

#
# if true, build with very strict error checking (disabled in src release tarballs)
#
set(DEVELOPER_MODE true)

#
# setup project/version
#

set(THIS_PROJECT_NAME "manta")

project (${THIS_PROJECT_NAME})

# find interpreters
find_package(PythonInterp)
if (NOT PYTHONINTERP_FOUND)
    message (WARNING "No python interpreter found, disabling optional python build and installation components. Installed workflow requires python interpreter to run")
endif()

# paths:
set(THIS_REDIST_DIR "${CMAKE_SOURCE_DIR}/../redist")
set(THIS_MODULE_DIR "${CMAKE_SOURCE_DIR}/cmake")
set(THIS_SOURCE_QC_DIR "${CMAKE_SOURCE_DIR}/srcqc")

# Create package versioning target - version derived from git describe except for
# release tarballs. Version file is treated as always out of date.
set(THIS_BUILDTIME_CONFIG_FILE "${CMAKE_CURRENT_BINARY_DIR}/buildTimeConfigInfo.txt")
set(THIS_BUILDTIME_CONFIG_TARGET "${THIS_PROJECT_NAME}_buildtime_config")
add_custom_target(
    ${THIS_BUILDTIME_CONFIG_TARGET}
    ALL
    COMMAND ${CMAKE_COMMAND}
    -D REDIST_DIR=${THIS_REDIST_DIR}
    -D CONFIG_FILE=${THIS_BUILDTIME_CONFIG_FILE}
    -P ${THIS_MODULE_DIR}/getBuildTimeConfigInfo.cmake
    COMMENT "Updating buildtime config info")

set (THIS_GLOBALS_CMAKE "${THIS_MODULE_DIR}/globals.cmake")
set (THIS_MACROS_CMAKE "${THIS_MODULE_DIR}/macros.cmake")

# add custom functions:
include ("${THIS_MACROS_CMAKE}")

#
# setup build type
#
set(CMAKE_BUILD_TYPES Debug Release RelWithDebInfo ASan)
join(CMAKE_BUILD_TYPES ", " CMAKE_BUILD_TYPES_STRING)

if (NOT CMAKE_BUILD_TYPE)
    set(DEFAULT_CMAKE_BUILD_TYPE Release)
    set(CMAKE_BUILD_TYPE ${DEFAULT_CMAKE_BUILD_TYPE} CACHE STRING
        "Choose the type of build, options are: {${CMAKE_BUILD_TYPES_STRING}} (default: ${DEFAULT_CMAKE_BUILD_TYPE})" FORCE)
endif()

# check that a valid build type has been selected:
list(FIND CMAKE_BUILD_TYPES ${CMAKE_BUILD_TYPE} TINDEX)
if(TINDEX EQUAL -1)
    message (FATAL_ERROR "Selected build type, '${CMAKE_BUILD_TYPE}', is not in the set of available types: {${CMAKE_BUILD_TYPES_STRING}}")
endif()

message (STATUS "BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

#
# setup cpack
#
if (CPACK_GENERATOR)
    message (STATUS "Configuring to produce the following package types: ${CPACK_GENERATOR}")
    set(CPACK_PACKAGE_VENDOR "Illumina")
    #set(CPACK_PACKAGE_VERSION "${WORKFLOW_VERSION}")
    set(CPAKC_PACKAGE_CONTACT "csaunders@illumina.com")

    if (CPACK_GENERATOR STREQUAL "RPM")
        set(CPACK_RPM_PACKAGE_RELOCATABLE TRUE)
        set(CPACK_RPM_PACKAGE_LICENSE "Illumina Open Source Software License 1")
        set(CPACK_RPM_PACKAGE_REQUIRES "python >= 2.4.0")
    else ()
        set(CPACK_SET_DESTDIR ON)
    endif ()

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

# paths for binary components -- we don't currently change this from the non-binary components
set(THIS_EXEC_PREFIX "${CMAKE_INSTALL_PREFIX}")
set(THIS_BINDIR "${THIS_EXEC_PREFIX}/bin")
set(THIS_LIBDIR "${THIS_EXEC_PREFIX}/lib")
set(THIS_LIBEXECDIR "${THIS_EXEC_PREFIX}/libexec")

# paths w/o architecture dependency
set(THIS_PYTHON_LIBDIR "${CMAKE_INSTALL_PREFIX}/lib/python")
set(THIS_DATADIR "${CMAKE_INSTALL_PREFIX}/share")
set(THIS_DOCDIR "${CMAKE_INSTALL_PREFIX}/doc")
set(THIS_DEMODIR "${THIS_DATADIR}/demo")

message (STATUS "install prefix: ${CMAKE_INSTALL_PREFIX}")

set(THIS_OPT "${THIS_PROJECT_NAME}_opt")
add_custom_target(${THIS_OPT} ALL)

# redist includes
include ("${THIS_MODULE_DIR}/boost.cmake")

if (NOT WIN32)
    # global source QC -- don't allow non-ascii chars in source files:
    set(THIS_SOURCE_CHECK "${THIS_PROJECT_NAME}_source_check")
    add_custom_target(${THIS_SOURCE_CHECK}
        ALL
        COMMAND bash ${THIS_SOURCE_QC_DIR}/check_for_nonascii_source.bash
        )

    # force this to run early:
    add_dependencies(${THIS_OPT} ${THIS_SOURCE_CHECK})
endif ()

# force unit tests at compile time, and set them to be more verbose than
# cmake's default "test" target:
set(THIS_UNITTESTS "${THIS_PROJECT_NAME}_unittests")
add_custom_target(${THIS_UNITTESTS}
    ALL
    COMMAND ${CMAKE_CTEST_COMMAND} --force-new-ctest-process --output-on-failure
    )

add_subdirectory (cmake/preInstall)
add_subdirectory (opt)
add_subdirectory (c++)
add_subdirectory (python)
if (NOT WIN32)
    add_subdirectory (demo)
    add_subdirectory (markdown)
endif ()
add_subdirectory (cmake/postInstall)
