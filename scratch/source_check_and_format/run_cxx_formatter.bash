#!/usr/bin/env bash
#
# run astyle on all cxx source
#

set -o nounset

scriptName=$(basename $0)

if ! which -a astyle > /dev/null 2>&1 ; then
    echo "ERROR: Can't find required utility 'astyle' in PATH" 1>&2
    exit 1;
fi

thisDir=$(dirname $0)

cxx_base_dir=$thisDir/../../src/c++


cd $cxx_base_dir 
astyle \
--style=ansi \
--align-pointer=type \
--max-instatement-indent=80 \
--min-conditional-indent=0 \
--pad-header \
--lineend=linux \
--suffix=none \
--recursive \
"*.cpp" "*.hh" "*.h" "*.h.in"

#--keep-one-line-blocks \
#--keep-one-line-statements \
