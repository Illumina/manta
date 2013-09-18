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
## CMake configuration file for all the c++ executables
##
## author Roman Petrovski
##
################################################################################

include (${MANTA_GLOBALS_CMAKE})

get_filename_component(CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
message (STATUS "Adding the c++    program subdirectory: ${CURRENT_DIR_NAME}")
include_directories (BEFORE SYSTEM ${MANTA_CXX_BEFORE_SYSTEM_INCLUDES})
include_directories (${MANTA_CXX_ALL_INCLUDES})
include_directories (${CMAKE_CURRENT_BINARY_DIR})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories (${MANTA_CXX_CONFIG_H_DIR})
