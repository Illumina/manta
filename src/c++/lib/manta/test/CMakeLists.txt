#
# Manta
# Copyright (c) 2013-2014 Illumina, Inc.
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

#include_directories (BEFORE SYSTEM "${SAMTOOLS_DIR}")
set(ADDITIONAL_UNITTEST_LIB manta_svgraph manta_blt_util)
include(${MANTA_CXX_TEST_LIBRARY_CMAKE})
