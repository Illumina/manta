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
## Configuration file for the unit tests subdirectory
##
## author Ole Schulz-Trieglaff
##
################################################################################

set(IS_QUIET true)
include(${THIS_CXX_EXECUTABLE_CMAKE})

configure_files("${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}" "*.cpp")

file(GLOB TEST_SOURCE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(TEST_SOURCE_FILE ${TEST_SOURCE_FILES})
    set(TEST_SOURCE ${TEST_SOURCE} "${CMAKE_CURRENT_BINARY_DIR}/${TEST_SOURCE_FILE}")
endforeach()

set(TEST_TARGET_NAME "${THIS_PROJECT_NAME}_unit_test_${THIS_LIB_DIR}")

add_executable(${TEST_TARGET_NAME} ${TEST_SOURCE})
add_dependencies(${TEST_TARGET_NAME} ${THIS_OPT})

if (THIS_LIBRARY_SOURCES)
    set(ADDITIONAL_UNITTEST_LIB ${ADDITIONAL_UNITTEST_LIB} ${THIS_PROJECT_NAME}_${THIS_LIB_DIR})
endif ()

target_link_libraries (${TEST_TARGET_NAME} ${ADDITIONAL_UNITTEST_LIB} ${THIS_AVAILABLE_LIBRARIES}
                      ${SAMTOOLS_LIBRARY} ${HTSLIB_LIBRARY} ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES} ${THIS_ADDITIONAL_LIB})

set(TEST_BINARY ${CMAKE_CURRENT_BINARY_DIR}/${TEST_TARGET_NAME})

add_test(${TEST_TARGET_NAME} ${TEST_BINARY} "--log_level=test_suite")
add_dependencies(${THIS_UNITTESTS} ${TEST_TARGET_NAME})
