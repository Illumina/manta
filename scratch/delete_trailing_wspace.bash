#!/usr/bin/env bash

set -o nounset
set -o pipefail

reltoabs() {
    (cd $1; pwd -P)
}

scriptdir=$(dirname $0)
basedir=$(reltoabs $scriptdir/..)


if [ $# != 1 ] || [ "$1" != "-imeanit" ]; then
    echo
    echo "usage: $0 -imeanit"
    echo ""
    echo "Cleanup whitespace in all manta code"
    echo "This is a slightly dangerous script, make sure your work is committed before running it"
    echo
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
        -name "*.py"
}

get_source() {
    for f in $basedir/src/*; do
        dir=$(basename $f)
        if [ $dir == "submodule" ]; then continue; fi
        find_cmake_source $f
        find_cxx_source $f
        find_script_source $f
    done
}

for f in $(get_source); do
    echo "checking: $f"
    sed -i 's/[ 	]*$//' $f
done
