#!/usr/bin/env bash

#
# facilitate comparison of manta output from two runs by stripping out the header and Manta id tags
#

set -o nounset
set -o pipefail

scriptName=$(basename $0)

#
# unzip or cat file as required:
#
optionalUngzip() {
    infile=$1
    if file -b $infile | grep -q gzip; then
        gzip -dc $infile
    else
        cat $infile
    fi
}

filterHeader() {
    awk '!/^#/'
}

stripMantaIds() {
    sed "s/Manta[^0-9]*:[0-9]*:[0-9]*:[0-9]*:[0-9]*:[0-9]*\(:[0-9]*\)\?//g"
}


#
# Optional add-on filters
#
stripPairCounts() {
    sed "s/PAIR_COUNT=[0-9]*//g"
}

stripQual() {
    awk 'BEGIN {FS="\t"; OFS="\t";} {$6=""; print}'
}

stripSample() {
    awk 'BEGIN {FS="\t"; OFS="\t";} {for(i=10;i<=NF;++i) { $i=""; } print;}'
}


#
# optionally ungzip, then remove header and remove IDs from each vcf
#
stripVcfCore() {
    infile=$1

    optionalUngzip $infile |\
    filterHeader |\
    stripMantaIds
}


#
# strip vcf and add additional info to make diff more informative:
#
stripVcf() {
    infile=$1

    # print input filename so that it's easy to figure out the diff polarity and see global lineCount diff:
    echo "$scriptName filteredVcf: $infile"
    lineCount=$(stripVcfCore $infile | wc -l)
    echo "$scriptName variantLineCount: $lineCount"

    # extra spaces keep the filename/lineCount diff above from attaching to a change on line 1 of the file:
    echo -e "\n\n\n"

    stripVcfCore $infile
}



#
# parse cmdline and diff:
#
if [ $# != 2 ]; then
    cat <<END
usage: $0 file1.vcf[.gz] file2.vcf[.gz]

This script helps to compare two manta vcf files. The vcfs can be
bgziped, gziped or uncompressed. Each file will be uncomressed and the
header + all ID Manta* fields will be stripped out.
END
    exit 2
fi

file1=$1
file2=$2

diff -U 0 <(stripVcf $file1) <(stripVcf $file2)
