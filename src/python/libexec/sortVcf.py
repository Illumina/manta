#!/usr/bin/env python
#
# Manta
# Copyright (c) 2013 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
#

"""
sort input vcf
"""

import os, sys
import re



def getKeyVal(string,key) :
    match=re.search("%s=([^;]*);?" % (key) ,string)
    if match is None : return None
    return match.group(1);


VCF_CHROM = 0
VCF_POS = 1
VCF_REF = 3
VCF_ALT = 4
VCF_QUAL = 5
VCF_FILTER = 6
VCF_INFO = 7



class VcfRecord :
    def __init__(self, line, isUnique) :
        self.line = line
        w=line.strip().split('\t')
        self.chrom=w[VCF_CHROM]
        self.pos=int(w[VCF_POS])
        if isUnique :
            self.ref=w[VCF_REF]
            self.alt=w[VCF_ALT]
            self.qual=w[VCF_QUAL]
            self.isPass=(w[VCF_FILTER] == "FILTER")
        self.endPos=self.pos+len(w[VCF_REF])-1
        val = getKeyVal(w[VCF_INFO],"END")
        if val is not None :
            self.endPos = int(val)



def processFile(isUnique, vcfFile, isFirst, header, recList) :
    """
    read in a vcf file
    """

    for line in open(vcfFile) :
        if line[0] == "#" :
            if isFirst : header.append(line)
        else :
            recList.append(VcfRecord(line, isUnique))



def getOptions() :

    from optparse import OptionParser

    usage = "usage: %prog [options] [vcf [vcf...]] > sorted_vcf"
    parser = OptionParser(usage=usage)

    parser.add_option("-u", dest="isUnique",action="store_true",default=False,
                      help="filter all but one record with the same {CHR,POS,REF,ALT}")

    (options,args) = parser.parse_args()

    if len(args) == 0 :
        parser.print_help()
        sys.exit(2)

    # validate input:
    for arg in args :
        if not os.path.isfile(arg) :
            raise Exception("Can't find input vcf file: " +arg)

    return (options,args)



def resolveRec(recEqualSet,recList) :
    """
    determine which of a set of 'equal' vcf records is the best
    right now best is a record with PASS in the filter field, and secondarily the high quality
    """

    if not recEqualSet: return

    bestIndex=0
    bestQual=0.
    bestIsPass=False
    for (index,rec) in enumerate(recEqualSet) :
        try:
            rec.qual = float(rec.qual)
        except ValueError:
            rec.qual = 0.

        assert rec.qual >= 0.

        isNewPass=((not bestIsPass) and rec.isPass)
        isHighQual=((bestIsPass == rec.isPass) and (rec.qual > bestQual))
        if (isNewPass or isHighQual) :
            bestIndex = index
            bestQual = rec.qual
            bestIsPass = rec.isPass

    recList.append(recEqualSet[bestIndex])



def main() :

    outfp = sys.stdout

    (options,args) = getOptions()

    header=[]
    recList=[]

    isFirst=True
    for arg in args :
        processFile(options.isUnique, arg, isFirst, header, recList)
        isFirst-False

    recList.sort(key = lambda x: (x.chrom, x.pos, x.endPos))

    for line in header :
        outfp.write(line)

    if options.isUnique :
        recList2 = []
        recEqualSet = []
        lastRec = None
        for vcfrec in recList :
            rec = (vcfrec.chrom,vcfrec.pos,vcfrec.ref,vcfrec.alt)
            if rec != lastRec :
                resolveRec(recEqualSet,recList2)
                recEqualSet = []
            recEqualSet.append(vcfrec)
            lastRec = rec
        resolveRec(recEqualSet,recList2)
        recList = recList2

    for vcfrec in recList :
        outfp.write(vcfrec.line)


main()
