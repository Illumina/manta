#
# Manta - Structural Variant and Indel Caller
# Copyright (c) 2013-2015 Illumina, Inc.
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
## file CMakeLists.txt
##
################################################################################

# run at build time
set(FA_FILE "Homo_sapiens_assembly19.COST16011_region.fa")
execute_process(
    COMMAND bzip2 -dc ${CMAKE_CURRENT_SOURCE_DIR}/${FA_FILE}.bz2
    OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/${FA_FILE}
    )

file (GLOB DATA_LIST *)
list (REMOVE_ITEM DATA_LIST ${CMAKE_CURRENT_SOURCE_DIR}/${FA_FILE}.bz2)
list (REMOVE_ITEM DATA_LIST ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt)
list (APPEND DATA_LIST ${CMAKE_CURRENT_BINARY_DIR}/${FA_FILE})

install(FILES ${DATA_LIST} DESTINATION ${THIS_MANTADEMODIR}/data)

