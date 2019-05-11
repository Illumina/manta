#!/usr/bin/env bash
#
# run clang-format on all cxx source
#

set -o nounset

this_dir=$(dirname $0)

cxx_base_dir=$this_dir/../../src/c++

get_c_and_cpp_files() {
  find $cxx_base_dir \( -name *.cpp -or -name *.cc -or -name *.c -or -name *.hpp -or -name *.hh -or -name *.h \)
}

# remove windows line endings from source:
get_c_and_cpp_files | xargs -P8 -n1 sed $'s/\r$//' -i

# general c/c++ source reformatting:
get_c_and_cpp_files | xargs -P8 -n1 $this_dir/clang-format -style=file -i
 
