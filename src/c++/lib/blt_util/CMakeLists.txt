#
# Manta
# Copyright (c) 2013 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/downloads/sequencing/licenses/>.
#

###############################################################################
##
## Configuration file for the c++/blt_util subfolder
##
## author Come Raczy
##
################################################################################

# samtools used by classes: bam_streamer and bam_record
#
include_directories (BEFORE SYSTEM "${SAMTOOLS_DIR}")
include(${MANTA_CXX_LIBRARY_CMAKE})
