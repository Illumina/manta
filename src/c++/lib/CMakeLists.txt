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
## Configuration file for the lib subfolder
##
## author Ole Schulz-Trieglaff
##
################################################################################

# Library build helper macro
macro(buildLibraries library_dirs library_targets)
    # Set target names
    set (${library_targets} "")
    foreach (library_dir ${library_dirs})
        get_library_target_name(${library_dir} LIBRARY_TARGET_NAME)
        list (APPEND ${library_targets} ${LIBRARY_TARGET_NAME})
    endforeach ()
    set (${library_targets} ${${library_targets}} PARENT_SCOPE)

    # Process library subdirectories
    foreach (library_dir ${library_dirs})
        add_subdirectory(${library_dir})
    endforeach ()
endmacro()


#
# Build libraries needed for unit testing separately
#
# This is logic required for unit-testing, but not deployed production code. This must be run before the
# primary loop below, so that PROJECT_TEST_LIBRARY_TARGETS is populated prior to configuring the unit test
# directories.
#
set (PROJECT_TEST_LIBRARY_DIRS test)
buildLibraries("${PROJECT_TEST_LIBRARY_DIRS}" PROJECT_TEST_LIBRARY_TARGETS)


#
# Build all primary libraries for the project
#
# These are libraries needed by most installed project applications.
#
# note that order of PROJECT_PRIMARY_LIBRARY_DIRS approximates dependency chain, where
# libraries are listed before libraries upon which they depend:
set (PROJECT_PRIMARY_LIBRARY_DIRS manta svgraph format alignment assembly options appstats htsapi common blt_util)
buildLibraries("${PROJECT_PRIMARY_LIBRARY_DIRS}" PROJECT_PRIMARY_LIBRARY_TARGETS)


#
# handle applications directory separately
#
add_subdirectory(applications)
