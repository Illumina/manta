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
## Configuration file for the lib subfolder
##
## author Ole Schulz-Trieglaff
##
################################################################################

##
## List of libraries
##
## note that order approximates dependency chain, where libraries are listed
## after libraries upon which they depend:
set (MANTA_ALL_LIBRARIES blt_util common options assembly alignment format svgraph manta)

##
## Build all the libraries for the project
## MANTA_AVAILABLE_LIBRARIES is incrementally updated
##
set (MANTA_AVAILABLE_LIBRARIES "")
foreach (MANTA_LIB_DIR ${MANTA_ALL_LIBRARIES})
    add_subdirectory(${MANTA_LIB_DIR})

    set(MANTA_AVAILABLE_LIBRARIES manta_${MANTA_LIB_DIR} ${MANTA_AVAILABLE_LIBRARIES})
endforeach ()

set (MANTA_AVAILABLE_LIBRARIES ${MANTA_AVAILABLE_LIBRARIES} PARENT_SCOPE)


#
# handle applications separately
#
add_subdirectory(applications)

