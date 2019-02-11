#
# Manta - Structural Variant and Indel Caller
# Copyright (c) 2013-2019 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

################################################################################
##
## Configuration file for the unit tests subdirectory
##
## author Ole Schulz-Trieglaff
##
################################################################################

# Incomplete list of required settings:
# LIBRARY_DIR_NAME directory name of parent library
# LIBRARY_TARGET_NAME cmake target name of parent library

include (${THIS_CXX_COMMMON_CMAKE})

setup_testConfig()

set(TEST_TARGET_NAME "${THIS_PROJECT_NAME}_unit_test_${LIBRARY_DIR_NAME}")

if (WIN32)
    # create a fake library target on win32 instead of linking and running the unit test
    # this creates a project in VS that allows for interaction with the unit test code

    # add all files to TEST_SOURCE to make the IDE project more usable:
    file(GLOB TMP_TARGET_FILES *)

    add_library     (${TEST_TARGET_NAME} STATIC ${TMP_TARGET_FILES})
    #set_source_files_properties(thefile PROPERTIES HEADER_FILE_ONLY TRUE)
    add_dependencies(${TEST_TARGET_NAME} ${THIS_OPT})
else ()
    file(GLOB TEST_SOURCE *.cpp)
    add_executable(${TEST_TARGET_NAME} ${TEST_SOURCE})
    add_dependencies(${TEST_TARGET_NAME} ${THIS_OPT})

    target_link_libraries (${TEST_TARGET_NAME} ${LIBRARY_TARGET_NAME}
                           ${PROJECT_TEST_LIBRARY_TARGETS} ${PROJECT_PRIMARY_LIBRARY_TARGETS}
                           ${HTSLIB_LIBRARY} ${Boost_LIBRARIES} ${THIS_ADDITIONAL_LIB})

    set(TEST_BINARY ${CMAKE_CURRENT_BINARY_DIR}/${TEST_TARGET_NAME})

    add_test(${TEST_TARGET_NAME} ${TEST_BINARY} "--log_level=test_suite")
endif ()

# make the target project use folders when applying cmake IDE generators like Visual Studio
file(RELATIVE_PATH THIS_RELATIVE_LIBDIR "${THIS_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET ${TEST_TARGET_NAME} PROPERTY FOLDER "${THIS_RELATIVE_LIBDIR}")

add_dependencies(${THIS_UNITTESTS} ${TEST_TARGET_NAME})
