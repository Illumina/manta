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
set (THIS_ALL_LIBRARIES blt_util common htsapi appstats options assembly alignment format svgraph truth manta)

##
## Build all the libraries for the project
## THIS_AVAILABLE_LIBRARIES is incrementally updated
##
set (THIS_AVAILABLE_LIBRARIES "")
foreach (THIS_LIB_DIR ${THIS_ALL_LIBRARIES})
    add_subdirectory(${THIS_LIB_DIR})

    set(THIS_AVAILABLE_LIBRARIES ${THIS_PROJECT_NAME}_${THIS_LIB_DIR} ${THIS_AVAILABLE_LIBRARIES})
endforeach ()

set (THIS_AVAILABLE_LIBRARIES ${THIS_AVAILABLE_LIBRARIES} PARENT_SCOPE)


#
# handle applications separately
#
add_subdirectory(applications)

