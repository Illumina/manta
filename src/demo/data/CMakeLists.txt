#
# Starka
# Copyright (c) 2009-2014 Illumina, Inc.
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
## file CMakeLists.txt
##
################################################################################

# run at cmake time
set(FA_FILE "chr20_860k_only.fa")
execute_process(
    COMMAND bzip2 -dc ${CMAKE_CURRENT_SOURCE_DIR}/${FA_FILE}.bz2
    OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/${FA_FILE}
    )

file (GLOB DATA_LIST *)
list (REMOVE_ITEM DATA_LIST ${CMAKE_CURRENT_SOURCE_DIR}/${FA_FILE}.bz2)
list (REMOVE_ITEM DATA_LIST ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt)
list (APPEND DATA_LIST ${CMAKE_CURRENT_BINARY_DIR}/${FA_FILE})

install(FILES ${DATA_LIST} DESTINATION ${THIS_STRELKADEMODIR}/data)

