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
## Configuration file for the unit tests subdirectory
##
## author Ole Schulz-Trieglaff
##
################################################################################

set(IS_QUIET true)
include(${MANTA_CXX_EXECUTABLE_CMAKE})

file (GLOB TEST_SOURCE "*.cpp")

set(TEST_TARGET unit_test_${MANTA_LIB_DIR})

add_executable(${TEST_TARGET} ${TEST_SOURCE})
add_dependencies(${TEST_TARGET} MANTA_OPT)

if (MANTA_LIBRARY_SOURCES)
    set(ADDITIONAL_UNITTEST_LIB ${ADDITIONAL_UNITTEST_LIB} manta_${MANTA_LIB_DIR})
endif ()

target_link_libraries (${TEST_TARGET} ${ADDITIONAL_UNITTEST_LIB} ${MANTA_AVAILABLE_LIBRARIES}
                      ${SAMTOOLS_LIBRARY} ${TABIX_LIBRARY} ${Boost_LIBRARIES} ${MANTA_ADDITIONAL_LIB})

set(TEST_BINARY ${CMAKE_CURRENT_BINARY_DIR}/${TEST_TARGET})

add_test(${TEST_TARGET} ${TEST_BINARY} "--log_level=test_suite")
add_dependencies(MANTA_UNITTESTS ${TEST_TARGET})
