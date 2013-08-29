#!/usr/bin/env bash


#
# improve comparison of manta output from two runs by stripping out the header and MantaBND id tags
#

set -o nounset

stripVcf() {
   awk '!/^#/' | sed "s/MantaBND:[0-9]*:[0-9]*:[0-9]*:[0-9]*:[0-9]*//g" 
}

file1=$1
file2=$2

diff <(stripVcf < $file1) <(stripVcf < $file2)

