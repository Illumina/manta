#!/usr/bin/env bash

#
# This documents the conversion of full boost to the subset contained in this project
#
# The subset used here only removes doc and support items to reduce the tarball size -- all
# code/libraries are present
#

set -o nounset
set -o xtrace

rel2abs() {
  pwd -P $1
}

thisDir=$(rel2abs $(dirname $0))

mkdir -p output
cd output

boost_name=boost_1_58_0
output_name=${boost_name}

boost_tarball=$thisDir/$boost_name.tar.bz2

if ! [ -f $boost_tarball ]; then
    echo "Can't find input boost tarball '$boost_tarball'"
    exit 1
fi

tar -xjf $boost_tarball

for ddir in doc more status; do
    rm -rf $output_name/$ddir
done

# remove docs:
(
cd $output_name
find . -name doc -type d -print | xargs rm -rf
)

# tarball up:

tar -c $output_name -f - | bzip2 -c -9 >| $output_name.tar.bz2
