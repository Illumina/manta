#!/usr/bin/env bash

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

