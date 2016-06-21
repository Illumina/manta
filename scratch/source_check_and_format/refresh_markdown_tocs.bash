#!/usr/bin/env bash

set -o nounset
set -o pipefail

reltoabs() {
    (cd $1; pwd -P)
}

scriptdir=$(dirname $0)
basedir=$(reltoabs $scriptdir/../..)

if ! [ -d $basedir/scratch/util ]; then
    echo "ERROR: [$0] Can't find source repo root directory"
    exit 1
fi

tocMaker=$basedir/scratch/util/tocMaker.py


if [ $# != 0 ]; then
    cat <<EOF

usage: $0

Refresh TOCs in any markdown file containing a "# Table of Contents" or similar.

EOF
    exit 2
fi


find_markdown_source() {
    base_dir=$1
    find $base_dir -type f \
        -name "*.md"
}


get_source() {
    find_markdown_source $basedir
}

for f in $(get_source); do
    echo "updating_toc: $f"
    newF=$f.$$.toc
    $tocMaker --depth=3 < $f > $newF
    if cmp -s $f $newF; then
        rm $newF
    else
        mv $f $f.$$.tocbak
        mv $newF $f
    fi
done
