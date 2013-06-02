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

foreach (MANTA_DIR ${MANTA_TEST_DIRS})
    message (STATUS "Testing access to ${MANTA_DIR}...")
    execute_process(COMMAND bash -c "mkdir -p ${MANTA_DIR}/test && rmdir ${MANTA_DIR}/test" RESULT_VARIABLE TMP_RESULT )
    if (TMP_RESULT)
        message (STATUS "ERROR: Directory is not writeable: ${MANTA_DIR}")
        message (STATUS "If you don't have administrator access to the "
                         "target installation location, please use --prefix "
                         "command-line option when configuring iSAAC. "
                         "Please use configure --help for all installer "
                         "command-line options details.")
        message (FATAL_ERROR "ERROR: installation cannot continue")
    else ()
        message (STATUS "Directory is writeable: ${MANTA_DIR}")
    endif ()
endforeach ()
