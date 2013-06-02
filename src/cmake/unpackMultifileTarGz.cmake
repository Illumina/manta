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
## file unpackMultifileTarGz.cmake
##
## extracts the set of files matching the *.tar.gz.* pattern assuming they
## concatenate into a valid tar.gz file if ordered lexicographically
## arguments:
##  SOURCE_ARCHIVE_FOLDER - Location of *.tar.gz.* files
##  TARGET_FOLDER         - Destination folder. Tar will -C ${TARGET_FOLDER}
##  TAR_EXCLUDE_PATTERN   - will add --exclude ${TAR_EXCLUDE_PATTERN} to the 
##                          tar cmd line if not empty
##
## author Roman Petrovski
##
################################################################################

file(GLOB ARCHIVE_VOLUMES ${SOURCE_ARCHIVE_FOLDER}/*.tar.gz.*)

IF(NOT EXISTS ${TARGET_FOLDER}/.unpacked)

    IF(NOT "${TAR_EXCLUDE_PATTERN}" STREQUAL "")
        set(TAR_ARGS --exclude "${TAR_EXCLUDE_PATTERN}")
    ENDIF()

    execute_process(COMMAND mkdir -p ${TARGET_FOLDER})
    execute_process(COMMAND echo "${ARCHIVE_VOLUMES}"
                    COMMAND sed "s/;/\\n/g"
                    COMMAND sort
                    COMMAND xargs -L 1 cat
                    COMMAND tar -xzv ${TAR_ARGS} -C ${TARGET_FOLDER}
                    RESULT_VARIABLE UNPACK_RESULT
                   )
    IF(NOT ${UNPACK_RESULT})
        execute_process(COMMAND touch ${TARGET_FOLDER}/.unpacked)
    ENDIF()
ENDIF(
