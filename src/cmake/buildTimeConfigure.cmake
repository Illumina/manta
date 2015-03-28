#
# Starka
# Copyright (c) 2009-2014 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
#

# \author Chris Saunders

# update file with various build-time properties
# requires CONFIG_FILE SOURCE_FILE DEST_FILE
#
# expected CONFIG_FILE format is:
# "key1\tvalue1\n"
# "key2\tvalue2\n" ...
#

file (READ ${CONFIG_FILE} CONFIG_LINES)
STRING(REPLACE "\n" ";" CONFIG_LINES "${CONFIG_LINES}")
foreach (CONFIG_LINE ${CONFIG_LINES})
    STRING(REPLACE "\t" ";" CONFIG_LIST "${CONFIG_LINE}")
    list (GET CONFIG_LIST 0 PAIR0)
    list (GET CONFIG_LIST 1 PAIR1)
    set(${PAIR0} "${PAIR1}")
endforeach ()
configure_file(${SOURCE_FILE} ${DEST_FILE} @ONLY)
