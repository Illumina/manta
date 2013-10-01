#!/usr/bin/env bash
#
# check c++ source for "#define DEBUG"... provides quick check of code for debug output prior to release
#

set -o nounset


thisDir=$(dirname $0)

cxx_base_dir=$thisDir/../../src/c++

find_cxx_source() {
    base_dir=$1
    find $base_dir -type f \
        -name "*.cpp" -or \
        -name "*.c" -or \
        -name "*.hh" -or \
        -name "*.h"
}


is_error=false
for f in $(find_cxx_source $cxx_base_dir); do
    grep --color='auto' -n -H "^#define DEBUG" $f

    if [ $? != 1 ]; then
        is_error=true
    fi
done

if $is_error; then
    echo -e "\nERROR: source defines DEBUG preprocessor flags\n" 1>&2
    exit 1
fi

