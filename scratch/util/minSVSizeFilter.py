#!/usr/bin/env python
#
# Manta - Structural Variant and Indel Caller
# Copyright (c) 2013-2017 Illumina, Inc.
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

import sys
import re


def getKeyVal(string,key) :
    match=re.search("%s=([^;\t]*);?" % (key) ,string)
    if match is None : return None
    return match.group(1);


class VCFID :
    CHROM = 0
    POS = 1
    REF = 3
    ALT = 4
    QUAL = 5
    FILTER = 6
    INFO = 7



def getOptions() :

    from optparse import OptionParser

    usage = "usage: %prog [options] < in.vcf > out.vcf"
    parser = OptionParser(usage=usage)

    parser.add_option("--minSize", dest="minSize", type="int",
                      help="minimum indel size, no default (required)")

    (opt,args) = parser.parse_args()

    if (opt.minSize is None) or (len(args) != 0) :
        parser.print_help()
        sys.exit(2)

    # validate input:
    if opt.minSize < 1:
        raise Exception("Invalid minSize value: %i" % (opt.minSize))

    return (opt,args)


def main() :

    infp=sys.stdin
    outfp=sys.stdout

    (opt,args) = getOptions()


    for line in infp :
        if line[0] == "#" :
            outfp.write(line)
            continue

        w=line.strip().split('\t')
        assert(len(w) > VCFID.INFO)
    
        svtype = getKeyVal(w[VCFID.INFO],"SVTYPE")
        if (svtype is not None) and (svtype == "BND") :
            outfp.write(line)
            continue

        if w[VCFID.ALT] == "<INS>" :
            outfp.write(line)
            continue

        svlen = getKeyVal(w[VCFID.INFO],"SVLEN")
        if svlen is None :
            ref=w[VCFID.REF]
            alt=w[VCFID.ALT]
            if alt.find("<") != -1 : continue
            svlen = max(len(ref),len(alt)) - 1

        if svlen is None : continue
        if abs(int(svlen)) < opt.minSize : continue

        outfp.write(line)


main()
