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
## Configuration file for the c++ subfolder
##
## author Come Raczy
##
################################################################################

set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set (THIS_CXX_EXECUTABLE_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxExecutable.cmake")
set (THIS_CXX_LIBRARY_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxLibrary.cmake")
set (THIS_CXX_CONFIGURE_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxConfigure.cmake")
set (THIS_CXX_TEST_LIBRARY_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxTestLibrary.cmake")

include ("${THIS_CXX_CONFIGURE_CMAKE}")

##
## refresh git label at build-time
##
add_dependencies(${THIS_OPT} ${CXX_BUILDTIME_CONFIG_TARGET})

##
## run cppcheck if any cxx source has been updated:
##
file(GLOB_RECURSE ALL_CXX_FILES *.cpp *.cc *.hh)

add_custom_command(
    OUTPUT cppcheck.done
    DEPENDS ${ALL_CXX_FILES}
    COMMAND ${PYTHON_EXECUTABLE} ${THIS_SOURCE_QC_DIR}/run_cppcheck.py ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Running c++ static analyzer")

if(${DEVELOPER_MODE})
    set(THIS_CPPCHECK "${THIS_PROJECT_NAME}_cppcheck")
    add_custom_target(${THIS_CPPCHECK} ALL DEPENDS cppcheck.done)

    add_dependencies(${THIS_OPT} ${THIS_CPPCHECK})
endif()

##
## Build all the libraries for the project
##
add_subdirectory (lib)


##
## build all the applications for the project
##
if (NOT WIN32)
    add_subdirectory (bin)
endif ()

##
## build the documentation when available
##
include  (FindDoxygen)
message (STATUS "Doxygen: ${DOXYGEN_EXECUTABLE}. Dot: ${DOXYGEN_DOT_EXECUTABLE}.")
if (DOXYGEN_FOUND)
    set (DOXYFILE ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    message (STATUS "Creating Doxygen config file: ${DOXYFILE}")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${DOXYFILE} @ONLY)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${DOXYFILE})
endif ()

