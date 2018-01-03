#!/usr/bin/env python
#
# Manta - Structural Variant and Indel Caller
# Copyright (c) 2013-2018 Illumina, Inc.
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
filter vcf to reduce inversions to a single vcf entry for cases where the same inversion is expressed by 
a left-right inversion combination.
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
        else :
            self.endPos = None
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

    usage = "usage: %prog < vcf > sorted_unique_inv_vcf"
    parser = OptionParser(usage=usage)

    (options,args) = parser.parse_args()

    if len(args) != 0 :
        parser.print_help()
        sys.exit(2)

    return (options,args)



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

    (options,args) = getOptions()

    header=[]
    recList=[]
    chromOrder=[]

    processStream(sys.stdin, chromOrder, header, recList)

    def vcfRecSortKey(x) :
        """
        sort vcf records for final output

        Fancy chromosome sort rules:
        if contig records are found in the vcf header, then sort chroms in that order
        for any chrom names not found in the header, sort them in lex order after the
        found chrom names
        """

        try :
            headerOrder = chromOrder.index(x.chrom)
        except ValueError :
            headerOrder = size(chromOrder)

        return (headerOrder, x.chrom, x.pos, x.endPos)

    recList.sort(key = vcfRecSortKey)

    for line in header :
        outfp.write(line)

    recList2 = []
    recEqualSet = []
    lastRec = None
    for vcfrec in recList :
        rec = (vcfrec.chrom, vcfrec.pos, vcfrec.endPos, vcfrec.svtype)
        if ((lastRec is None) or
          (rec[0] != lastRec[0]) or
          (rec[3] != "INV") or
          (lastRec[3] != "INV") or
          (abs(rec[1]-lastRec[1]) > 250) or
          (abs(rec[2]-lastRec[2]) > 250)) :
            resolveRec(recEqualSet,recList2)
            recEqualSet = []
   
        recEqualSet.append(vcfrec)
        lastRec = rec
    resolveRec(recEqualSet,recList2)
    recList = recList2

    for vcfrec in recList :
        outfp.write(vcfrec.line)


main()

