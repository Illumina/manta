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
## Configuration file for boost installation
##
## author Come Raczy
##
################################################################################

# simple helper for resetFindBoost
function(unsetall name)
    unset (${name} CACHE)
    unset (${name})
endfunction()


macro (resetFindBoost)

    set(BOOST_RESET_SYMBOLS FOUND INCLUDE_DIRS INCLUDE_DIR LIBRARIES LIBRARY_DIRS VERSION LIB_VERSION MAJOR_VERSION MINOR_VERSION SUBMINOR_VERSION USE_STATIC_LIBS USE_MULTITHREADED)

    foreach (TAG ${BOOST_RESET_SYMBOLS})
        unsetall (Boost_${TAG})
    endforeach()

    unset (ENV{BOOST_LIBRARYDIR})

    foreach (COMPONENT ${MANTA_BOOST_COMPONENTS})
        STRING(TOUPPER ${COMPONENT} UPPERCOMPONENT)
        unsetall (Boost_${UPPERCOMPONENT}_FOUND)
        unsetall (Boost_${UPPERCOMPONENT}_LIBRARY)
        unsetall (Boost_${UPPERCOMPONENT}_LIBRARY_RELEASE)
        unsetall (Boost_${UPPERCOMPONENT}_LIBRARY_DEBUG)
    endforeach ()

endmacro ()

#
# finds boost and sets variables so that
# it is being used for include and linking
#
macro(static_find_boost boost_version boost_components)

    # MANTA uses the static boost build... assuming pthread is not required in this case, but keeping check code around in case:
    #
    # pthread library required by boost
#    static_find_library(PTHREAD "pthread.h" pthread)
#    if    (HAVE_PTHREAD)
#        set  (MANTA_ADDITIONAL_LIB ${MANTA_ADDITIONAL_LIB} pthread)
#        message(STATUS "pthread supported")
#    else  ()
#        message(STATUS "pthread headers: ${PTHREAD_INCLUDE_DIR}")
#        message(STATUS "pthread library: ${PTHREAD_LIBRARY}")
#        message(FATAL_ERROR "pthread library is required to build")
#    endif ()

    find_package(Boost ${boost_version} REQUIRED ${boost_components})

    foreach(COMPONENT ${MANTA_BOOST_COMPONENTS})
        STRING(TOUPPER ${COMPONENT} UPPERCOMPONENT)
        set(HAVE_LIBBOOST_${UPPERCOMPONENT}       ${Boost_${UPPERCOMPONENT}_FOUND})
    endforeach()
endmacro()


if (MANTA_FORCE_STATIC_LINK)
    set(Boost_USE_STATIC_LIBS ON)
endif ()

find_package(Boost ${MANTA_BOOST_VERSION} COMPONENTS ${MANTA_BOOST_COMPONENTS})


# CMAKE_PARALLEL is only used if boost is found, but moving the setting here (outside of the if below) supresses a cmake warning:
if (NOT CMAKE_PARALLEL)
    set (CMAKE_PARALLEL "1")
endif ()

#
# If the right version of boost is not found, it will be built from the distribution
#
if (NOT Boost_FOUND)

    foreach(COMPONENT ${MANTA_BOOST_COMPONENTS})
        STRING(TOUPPER ${COMPONENT} UPPERCOMPONENT)
        if (${Boost_${UPPERCOMPONENT}_FOUND})
            set(FOUND_STATUS "found")
        else()
            set(FOUND_STATUS "NOT FOUND")
        endif()
        message(STATUS "Boost component: ${COMPONENT}\tstatus: ${FOUND_STATUS}")
    endforeach()

    if (BOOST_ROOT)
        message (STATUS "BOOST_ROOT is set to ${BOOST_ROOT} but boost ${MANTA_BOOST_VERSION} or one of its components was not found.")
        message (FATAL_ERROR "Unset BOOST_ROOT or set it to the root location of boost ${MANTA_BOOST_VERSION}.")
    endif()

    # Try to find it in target installation location
    resetFindBoost()
    message(STATUS "Boost ${MANTA_BOOST_VERSION} not found. Boost will be built from the distribution...")

    set(ENV{MANTA_BOOST_BUILD_COMPONENTS} "${MANTA_BOOST_BUILD_COMPONENTS}")
    set(ENV{MANTA_BOOST_VERSION} "${MANTA_BOOST_VERSION}")

    set(BOOST_BUILD_COMMAND bash "${MANTA_BOOTSTRAP_DIR}/installBoost.bash" "${MANTA_REDIST_DIR}" "${MANTA_BOOTSTRAP_DIR" "${CMAKE_PARALLEL}")

    string(REPLACE ";" " " PRETTY_BOOST_BUILD_COMMAND "${BOOST_BUILD_COMMAND}")
    message(STATUS "${PRETTY_BOOST_BUILD_COMMAND}")
    execute_process(COMMAND ${BOOST_BUILD_COMMAND} RESULT_VARIABLE TMP_RESULT)

    if (TMP_RESULT)
        message (FATAL_ERROR "Failed to build boost ${MANTA_BOOST_VERSION}")
    else ()
        message(STATUS "Successfuly built boost ${MANTA_BOOST_VERSION} from the distribution package...")
    endif ()

    set (BOOST_ROOT "${MANTA_BOOTSTRAP_DIR}")
endif ()

