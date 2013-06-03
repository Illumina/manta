#!/usr/bin/env bash
#
# if asytle is found, then run it on all cxx source
#

set -o nounset


if ! which -a astyle > /dev/null 2>&1 ; then exit 0; fi

thisDir=$(dirname $0)

cxx_base_dir=$thisDir/../src/c++


cd $cxx_base_dir 
astyle \
--align-pointer=type \
--max-instatement-indent=80 \
--min-conditional-indent=0 \
--pad-header \
--recursive \
*.cpp *.hh *.h

#--keep-one-line-blocks \
#--keep-one-line-statements \
