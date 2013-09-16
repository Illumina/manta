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
## author Roman Petrovski
##
################################################################################


message (STATUS "Verifying target directories access")

include ("${MANTA_GLOBALS_CMAKE}")
install(
    CODE "set(MANTA_TEST_DIRS \"\${MANTA_FULL_ETCDIR}\" \"\${MANTA_FULL_DATADIR}\"
                               \"\${MANTA_FULL_BINDIR}\" \"\${MANTA_FULL_LIBDIR}\"
                               \"\${MANTA_FULL_LIBEXECDIR}\")"
    SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/checkTargetPathWritable.cmake"
)
