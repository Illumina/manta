#!/usr/bin/env bash

#
# this script is used for changing license headers in the project, it attempts to cover all files
#

set -o nounset
set -o pipefail


rel2abs() {
    (cd $1; pwd -P)
}


thisdir=$(rel2abs $(dirname $0))


find_cxx_source() {
    base_dir=$1
    find $base_dir -type f \
        -name "*.cpp" -or \
        -name "*.hh"
}


find_python_source() {
    base_dir=$1
    find $base_dir -type f \
        -name "*.py"
}


find_cmake_source() {
    base_dir=$1
    find $base_dir -type f \
        -name "*.cmake" -or \
        -name "CMakeLists.txt"
}


find_shell_source() {
    base_dir=$1
    find $base_dir -type f \
        -name "*.bash" -or \
        -name "*.sh" -or \
        -name "configure"
}



reheader_file() {
    script="$1"
    file=$2

    if [ ! -f $file ]; then return; fi
    echo $file

    is_exe=false
    if [ -x $file ]; then
        is_exe=true
    fi
    tmpfile=$(mktemp)
    $script $thisdir/new_header < $file >| $tmpfile
    if [ $? != 0 ]; then echo "error on file $file"; exit 1; fi
    mv $tmpfile $file
    if [ $? != 0 ]; then echo "error on file $file"; exit 1; fi
    if $is_exe; then
        chmod +x $file
    fi
}


project_base_dir=$(rel2abs $thisdir/../../..)
cxx_base_dir=$project_base_dir/src/c++
python_base_dir=$project_base_dir/src
cmake_base_dir=$project_base_dir/src
shell_base_dir=$project_base_dir/src

for file in $(find_cxx_source $cxx_base_dir); do
    reheader_file "python $thisdir/reheader_cxx_file.py" $file
done

for file in $(find_python_source $python_base_dir) $(find_cmake_source $cmake_base_dir) $(find_shell_source $shell_base_dir); do
    reheader_file "python $thisdir/reheader_script_file.py" $file
done


