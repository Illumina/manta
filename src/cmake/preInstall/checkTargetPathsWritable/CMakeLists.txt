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
## author Roman Petrovski
##
################################################################################


message (STATUS "Verifying target directories access")

include ("${THIS_GLOBALS_CMAKE}")
install(
    CODE "set(THIS_TEST_DIRS \"\${THIS_FULL_DATADIR}\"
                               \"\${THIS_FULL_BINDIR}\" \"\${THIS_FULL_LIBDIR}\"
                               \"\${THIS_FULL_LIBEXECDIR}\")"
    SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/checkTargetPathWritable.cmake"
)
