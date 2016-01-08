#!/usr/bin/env bash
#
# Manta - Structural Variant and Indel Caller
# Copyright (c) 2013-2016 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

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
    sed "s/Manta[^0-9]*\(:[0-9]*\)\{5,7\}//g"
}


#
# Optional add-on filters
#
stripPairCounts() {
    sed "s/PAIR_COUNT=[0-9]*//g"
}

stripCigar() {
    sed "s/CIGAR=.*;\?//g"
}

stripQual() {
    awk 'BEGIN {FS="\t"; OFS="\t";} {$6=""; print}' | sed "s/JUNCTION_QUAL=[0-9]*//g" 
}

stripSample() {
    awk 'BEGIN {FS="\t"; OFS="\t";} {for(i=10;i<=NF;++i) { $i=""; } print;}'
}

stripSomScore() {
    sed "s/SOMATICSCORE=[0-9]*//g"
}

stripInvTags() {
    sed "s/INV5;//" | sed "s/INV3;//" | sed "s/INV5=1;//" | sed "s/INV3=1;//"
}

stripEvent() {
    sed "s/;EVENT=//"
}



#
# optionally ungzip, then remove header and remove IDs from each vcf
#
stripVcfCore() {
    infile=$1

    optionalUngzip $infile |\
    filterHeader |\
    stripMantaIds | stripSomScore
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
bgziped, gziped or uncompressed. Each file will be uncompressed and the
header + all ID Manta* fields will be stripped out.
END
    exit 2
fi

file1=$1
file2=$2

diff -U 0 <(stripVcf $file1) <(stripVcf $file2)
