#!/usr/bin/env bash

#
# This documents the conversion of full boost to the subset contained in this project
#
# The subset used here only removes doc and support items to reduce the tarball size -- all
# code/libraries are present
#

set -o nounset


boost_name=boost_1_49_0
subset_name=${boost_name}_subset


tar -xjf $boost_name.tar.bz2
mv $boost_name $subset_name

for ddir in doc more status; do
    rm -rf $subset_name/$ddir
done

# remove docs:
(
cd $subset_name
find . -name doc -type d -print | xargs rm -rf
)

# tarball up:

tar -c $subset_name -f - | bzip2 -c -9 >| $subset_name.tar.bz2
