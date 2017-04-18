#!/usr/bin/env bash
#
# if astyle is found, then run it on all cxx source
#

set -o nounset


if ! which -a astyle > /dev/null 2>&1 ; then exit 0; fi

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
