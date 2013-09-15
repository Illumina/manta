#!/usr/bin/env bash
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

# if cppcheck is found, run it on all c++ and return an error for *any* warning message:

set -o nounset
set -o pipefail

if ! which -a cppcheck > /dev/null 2>&1 ; then exit 0; fi

thisDir=$(dirname $0)

outFile=cppcheck.log


cppcheck \
--enable=all --std=c++03 --force --verbose --quiet \
--template='{file}:{line}:{severity}:{message}' \
--suppress=uninitMemberVar \
--suppress=unsignedLessThanZero \
--suppress=obsoleteFunctionsasctime \
--suppress=unusedFunction \
--suppress=unmatchedSuppression \
--suppress=missingInclude \
$thisDir 2>| $outFile

# xml output is usful for getting a warnings id field, which is what you need to supress it:
# --xml \

# this is more aggresive and includes more FPs
# --inconclusive \

filtered_output() {
    cat $outFile | grep -v "Cppcheck cannot find all the include files"
}

filtered_output 1>&2

# in cppcheck 1.59 missingInclude supression appears to be broken -- this is a workaround:
#
error_line_count=$(filtered_output | wc -l)

if [ "$error_line_count" != "0" ]; then exit 1; fi

touch cppcheck.done
