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
## Configuration file for boost installation
##
## author Come Raczy
##
################################################################################

# set to TRUE to see more detailed information about the boost find/build procedure:
set (DEBUG_FINDBOOST FALSE)
if (${DEBUG_FINDBOOST})
    set (Boost_DEBUG "ON")
    set (Boost_DETAILED_FAILURE_MSG "ON")
endif ()

macro (initBoostParams)
    # required boost libraries
    set (THIS_BOOST_VERSION 1.53.0)
    # note we default to alphabetical order here, except where boost libraries depend on
    # each other (timer->chrono)
    set (THIS_BOOST_COMPONENTS date_time filesystem program_options
                                regex serialization system timer chrono unit_test_framework)

    # the name given to boost.build and the library name are the same for all libraries, except
    # for test, so we need two lists now:
    set (THIS_BOOST_BUILD_COMPONENTS date_time filesystem program_options
                                     regex serialization system timer chrono test)
    set (Boost_USE_MULTITHREADED OFF)
    set (Boost_USE_STATIC_LIBS ON)
endmacro()

# simple helper for resetFindBoost
macro(unsetall name)
    unset (${name} CACHE)
    unset (${name})
endmacro()


macro (resetFindBoost)

    set(BOOST_RESET_SYMBOLS FOUND INCLUDE_DIRS INCLUDE_DIR LIBRARIES LIBRARY_DIRS VERSION LIB_VERSION MAJOR_VERSION MINOR_VERSION SUBMINOR_VERSION USE_STATIC_LIBS USE_MULTITHREADED)

    foreach (TAG ${BOOST_RESET_SYMBOLS})
        unsetall (Boost_${TAG})
    endforeach()

    unset (ENV{BOOST_LIBRARYDIR})

    foreach (COMPONENT ${THIS_BOOST_COMPONENTS})
        STRING(TOUPPER ${COMPONENT} UPPERCOMPONENT)
        unsetall (Boost_${UPPERCOMPONENT}_FOUND)
        unsetall (Boost_${UPPERCOMPONENT}_LIBRARY)
        unsetall (Boost_${UPPERCOMPONENT}_LIBRARY_RELEASE)
        unsetall (Boost_${UPPERCOMPONENT}_LIBRARY_DEBUG)
    endforeach ()

    initBoostParams()
endmacro ()


initBoostParams()

if (THIS_FORCE_STATIC_LINK)
    set(Boost_USE_STATIC_LIBS ON)
endif ()

set(BOOST_BOOTSTRAP_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/bootstrap/boost)
if (EXISTS "${BOOST_BOOTSTRAP_INSTALL_DIR}/boost_install_complete")
    set (BOOST_ROOT "${BOOST_BOOTSTRAP_INSTALL_DIR}")
endif ()

find_package(Boost ${THIS_BOOST_VERSION} COMPONENTS ${THIS_BOOST_COMPONENTS})

# CMAKE_PARALLEL is only used if boost is found, but moving the setting here (outside of the if below) supresses a cmake warning:
if (NOT CMAKE_PARALLEL)
    set (CMAKE_PARALLEL "1")
endif ()

#
# If the right version of boost is not found, it will be built from the distribution
#
if (NOT Boost_FOUND)
    foreach(COMPONENT ${THIS_BOOST_COMPONENTS})
        STRING(TOUPPER ${COMPONENT} UPPERCOMPONENT)
        if (${Boost_${UPPERCOMPONENT}_FOUND})
            set(FOUND_STATUS "found")
        else()
            set(FOUND_STATUS "NOT FOUND")
        endif()
        message(STATUS "Boost component: ${COMPONENT}\tstatus: ${FOUND_STATUS}")
    endforeach()

    if (BOOST_ROOT)
        message (STATUS "BOOST_ROOT is set to ${BOOST_ROOT} but boost ${THIS_BOOST_VERSION} or one of its components was not found.")
        message (FATAL_ERROR "Unset BOOST_ROOT or set it to the root location of boost ${THIS_BOOST_VERSION}.")
    endif()

    # Try to find it in target installation location
    resetFindBoost()
    message(STATUS "Boost ${THIS_BOOST_VERSION} not found. Boost will be built from the distribution...")

    set(ENV{THIS_BOOST_BUILD_COMPONENTS} "${THIS_BOOST_BUILD_COMPONENTS}")
    set(ENV{THIS_BOOST_VERSION} "${THIS_BOOST_VERSION}")

    set(THIS_BOOTSTRAP_DIR "${THIS_MODULE_DIR}/bootstrap")
    set(BOOST_BUILD_COMMAND bash "${THIS_BOOTSTRAP_DIR}/installBoost.bash" "${THIS_REDIST_DIR}" "${BOOST_BOOTSTRAP_INSTALL_DIR}" "${CMAKE_PARALLEL}")

    string(REPLACE ";" " " PRETTY_BOOST_BUILD_COMMAND "${BOOST_BUILD_COMMAND}")
    message(STATUS "${PRETTY_BOOST_BUILD_COMMAND}")
    execute_process(COMMAND ${BOOST_BUILD_COMMAND} RESULT_VARIABLE TMP_RESULT)

    if (TMP_RESULT)
        message (FATAL_ERROR "Failed to build boost ${THIS_BOOST_VERSION}")
    else ()
        message(STATUS "Successfuly built boost ${THIS_BOOST_VERSION} from the distribution package...")
    endif ()

    set (BOOST_ROOT "${BOOST_BOOTSTRAP_INSTALL_DIR}")
    find_package(Boost ${THIS_BOOST_VERSION} COMPONENTS ${THIS_BOOST_COMPONENTS})
endif ()

foreach(COMPONENT ${THIS_BOOST_COMPONENTS})
    STRING(TOUPPER ${COMPONENT} UPPERCOMPONENT)
    set(HAVE_LIBBOOST_${UPPERCOMPONENT} ${Boost_${UPPERCOMPONENT}_FOUND})
endforeach()
