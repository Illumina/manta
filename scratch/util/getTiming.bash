#!/usr/bin/env bash
#
# Manta - Structural Variant and Indel Caller
# Copyright (c) 2013-2017 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

#
# create a simple breakdown of task times for a manta run
#

set -o nounset
set -o pipefail

scriptName=$(basename $0)


#
# parse cmdline and diff:
#
if [ $# != 1 ]; then
    cat <<END
usage: $scriptName mantaRunDir

END
    exit 2
fi

dir=$1

if ! [ -d $dir ]; then
    echo "runDir $dir not found"
    exit 2
fi

elog=$dir/workspace/pyflow.data/logs/pyflow_tasks_stderr_log.txt

if ! [ -f $elog ]; then
    echo "runDir $dir does not contain expected files"
    exit 2
fi

cat $elog | awk '
/elapsedSec/ {
    if(/makeLocusGraph/) { graphSum += $9}
    else if(/generateCandidateSV/) { svSum += $9 }
    else { otherSum += $9 }
}
END {
    print "core-hours:";
    print "graph:",graphSum/(60*60);
    print "svcand:",svSum/(60*60);
    print "other:",otherSum/(60*60);
}'

