#!/usr/bin/env python
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

"""
remove intrachomosomal events above a specified size
"""

import os, sys
import re



def isInfoFlag(infoString,key) :
    word=infoString.split(";")
    for w in word :
        if w == key : return True
    return False


def getKeyVal(infoString,key) :
    match=re.search("%s=([^;\t]*);?" % (key) ,infoString)
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


class VcfRecord :
    def __init__(self, line) :
        self.line = line
        w=line.strip().split('\t')
        self.chrom=w[VCFID.CHROM]
        self.pos=int(w[VCFID.POS])
        self.qual=w[VCFID.QUAL]
        self.isPass=(w[VCFID.FILTER] == "PASS")
        self.Filter=w[VCFID.FILTER]
        self.endPos=self.pos+len(w[VCFID.REF])-1
        val = getKeyVal(w[VCFID.INFO],"END")
        if val is not None :
            self.endPos = int(val)
        val = getKeyVal(w[VCFID.INFO],"SOMATICSCORE")
        if val is not None :
            self.ss = int(val)
        else :
            self.ss = None
        self.svtype = getKeyVal(w[VCFID.INFO],"SVTYPE")
        self.isInv3 = isInfoFlag(w[VCFID.INFO],"INV3")
        self.isInv5 = isInfoFlag(w[VCFID.INFO],"INV5")



class Constants :

    import re

    contigpat = re.compile("^##contig=<ID=([^,>]*)[,>]")



def processStream(vcfFp, chromOrder, header, recList) :
    """
    read in a vcf stream
    """

    import re

    for line in vcfFp :
        if line[0] == "#" :
            header.append(line)
            match = re.match(Constants.contigpat,line)
            if match is not None :
                chromOrder.append(match.group(1))
        else :
            recList.append(VcfRecord(line))



def getOptions() :

    from optparse import OptionParser

    usage = "usage: %prog [options] < vcf > filtered_vcf"
    parser = OptionParser(usage=usage)
    
    parser.add_option("--maxSize", type="int",dest="maxSize",
                      help="maximum intrachrom event size [required] (no default)")

    (opt,args) = parser.parse_args()

    if len(args) != 0 :
        parser.print_help()
        sys.exit(2)

    if opt.maxSize is None :
        parser.print_help()
        sys.exit(2)

    return (opt,args)



def resolveRec(recEqualSet, recList) :
    """
    determine which of a set of vcf records presumed to refer to the same inversion are kept
    right now best is a record with PASS in the filter field, and secondarily the high quality
    """

    if not recEqualSet: return

    bestIndex=0
    bestSS=0.
    bestPos=0
    bestIsPass=False
    for (index,rec) in enumerate(recEqualSet) :
        assert rec.pos > 0

        isNewPass=((not bestIsPass) and rec.isPass)
        isHighQual=((bestIsPass == rec.isPass) and (rec.pos < bestPos)) #(rec.ss > bestSS))
        if (isNewPass or isHighQual) :
            bestIndex = index
            bestPos = rec.pos
            bestIsPass = rec.isPass

# potentially could reward two non-pass inversion calls here:
#    if not bestIsPass and (len(recEqualSet) == 2) :
#        if (recEqualSet[0].isInv3 and reEqualSet[1].isInv5) or
#            recEqualSet[1].isInv3 and reEqualSet[0].isInv5)) :

    recList.append(recEqualSet[bestIndex])



def main() :

    outfp = sys.stdout

    (opt,args) = getOptions()

    header=[]
    recList=[]
    chromOrder=[]

    processStream(sys.stdin, chromOrder, header, recList)

    for line in header :
        outfp.write(line)

    for vcfrec in recList :
        if (vcfrec.endPos-vcfrec.pos) > opt.maxSize : continue
        outfp.write(vcfrec.line)



main()

