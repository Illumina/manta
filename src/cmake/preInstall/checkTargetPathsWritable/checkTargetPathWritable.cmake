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

foreach (MANTA_DEST_DIR ${MANTA_DEST_DIRS})
    message (STATUS "Testing access to ${MANTA_DEST_DIR}...")
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "make_directory" "${MANTA_DEST_DIR}/test" RESULT_VARIABLE TMP_RESULT )
    if (TMP_RESULT)
        message (STATUS "ERROR: Directory is not writeable: ${MANTA_DEST_DIR}")
        message (STATUS "If you don't have administrator access to the "
                         "target installation location, please use --prefix "
                         "command-line option when configuring iSAAC. "
                         "Please use configure --help for all installer "
                         "command-line options details.")
        message (FATAL_ERROR "ERROR: iSAAC installation cannot continue")
    else ()
        message (STATUS "Directory is writeable: ${MANTA_DEST_DIR}")
    endif ()
endforeach ()
