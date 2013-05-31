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

################################################################################
##
## author Roman Petrovski
##
################################################################################


message (STATUS "Verifying target directories access")

include ("${MANTA_GLOBALS_CMAKE}")
install(
    CODE "set(MANTA_DEST_DIRS \"\${MANTA_DEST_ETCDIR}\" \"\${MANTA_DEST_DATADIR}\" 
                               \"\${MANTA_DEST_BINDIR}\" \"\${MANTA_DEST_LIBDIR}\" 
                               \"\${MANTA_DEST_LIBEXECDIR}\")"
    SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/checkTargetPathWritable.cmake"
)
