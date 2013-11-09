#!/usr/bin/env bash

#
# improve comparison of manta output from two runs by stripping out the header and Manta id tags
#

set -o nounset
set -o pipefail

scriptName=$(basename $0)

#
# unzip or cat file as required:
#
optionalUngzip() {
    infile=$1
    echo -e "$scriptName filtered vcf: $1\n\n\n"
    if file -b $infile | grep -q gzip; then
        gzip -dc $infile
    else
        cat $infile
    fi
}


#
# unzip, remove header and remove IDs from each vcf
#
stripVcf() {
   infile=$1

   optionalUngzip $1 |\
   awk '!/^#/' |\
   sed "s/Manta.*:[0-9]*:[0-9]*:[0-9]*:[0-9]*:[0-9]*//g" 
}


file1=$1
file2=$2

diff -U 0 <(stripVcf $file1) <(stripVcf $file2)

