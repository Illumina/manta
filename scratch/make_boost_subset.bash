#!/usr/bin/env bash

#
# This documents the conversion of full boost to the subset contained in this project
#
# The subset used here only removes doc and support items to reduce the tarball size -- all
# code/libraries are present
#

set -o nounset
set -o xtrace

mkdir -p output

boost_name=boost_1_56_0


tar -xjf $boost_name.tar.bz2
mv $boost_name output 
cd output

for ddir in doc more status; do
    rm -rf $boost_name/$ddir
done

# remove docs:
(
cd $boost_name
find . -name doc -type d -print | xargs rm -rf
)

# tarball up:

tar -c $boost_name -f - | bzip2 -c -9 >| $boost_name.tar.bz2
