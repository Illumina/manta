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

foreach (THIS_DIR ${THIS_TEST_DIRS})
    message (STATUS "Testing access to ${THIS_DIR}...")
    execute_process(COMMAND bash -c "mkdir -p ${THIS_DIR}/test && rmdir ${THIS_DIR}/test" RESULT_VARIABLE TMP_RESULT )
    if (TMP_RESULT)
        message (STATUS "ERROR: Directory is not writeable: ${THIS_DIR}")
        message (STATUS "If you don't have administrator access to the "
                         "target installation location, please use --prefix "
                         "command-line option during configuration. "
                         "Please see 'configure --help' for all installer "
                         "command-line options.")
        message (FATAL_ERROR "ERROR: installation cannot continue")
    else ()
        message (STATUS "Directory is writeable: ${THIS_DIR}")
    endif ()
endforeach ()
