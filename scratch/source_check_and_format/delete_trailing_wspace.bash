#!/usr/bin/env bash

set -o nounset
set -o pipefail

reltoabs() {
    (cd $1; pwd -P)
}

scriptdir=$(dirname $0)
basedir=$(reltoabs $scriptdir/../..)


if [ $# != 1 ] || [ "$1" != "-imeanit" ]; then
    cat <<EOF

usage: $0 -imeanit

Cleanup whitespace in all project code and ensure all files end in a newline.
This is a slightly dangerous script, make sure your work is committed before running it

EOF
    exit 2
fi


find_cmake_source() {
    base_dir=$1
    find $base_dir -type f \
        -name "*.cmake" -or \
        -name "CMakeLists.txt"
}

find_cxx_source() {
    base_dir=$1
    find $base_dir -type f \
        -name "*.cpp" -or \
        -name "*.c" -or \
        -name "*.hh" -or \
        -name "*.h"
}

find_script_source() {
    base_dir=$1
    find $base_dir -type f \
        -name "*.bash" -or \
        -name "*.sh" -or \
        -name "configure" -or \
        -name "*.py"
}

find_doc_source() {
    base_dir=$1
    find $base_dir -type f \
        -name "*.md" -or \
        -name "*.tex"
}


get_source() {
    for f in $basedir/src/*; do
        dir=$(basename $f)
        if [ $dir == "submodule" ]; then continue; fi
        find_cmake_source $f
        find_cxx_source $f
        find_script_source $f
    done
    find_cmake_source $basedir/redist
    find_doc_source $basedir
}

for f in $(get_source); do
    echo "checking: $f"
    sed -i 's/[ 	]*$//' $f
    python $scriptdir/ensureFileEndsInNewline.py $f
done
