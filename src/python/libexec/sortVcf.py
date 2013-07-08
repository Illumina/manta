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
# <https://github.com/downloads/sequencing/licenses/>.
#

"""
sort input vcf
"""

import os, sys


def getOptions() :

    from optparse import OptionParser

    usage = "usage: %prog [vcf [vcf...]] > sorted_vcf"
    parser = OptionParser(usage=usage)


    (options,args) = parser.parse_args()

    if len(args) == 0 :
        parser.print_help()
        sys.exit(2)

    # validate input:
    for arg in args :
        if not os.path.isfile(arg) :
            raise Exception("Can't find input vcf file: " +arg)

    return (options,args)



class VcfRecord :
    def __init__(self,line) :
        self.line = line
        w=line.strip().split('\t')
        self.chrom=w[0]
        self.pos=int(w[1])



def processFile(arg,isFirst,header,recList) :
    """
    read in a vcf file
    """

    for line in open(arg) :
        if line[0] == "#" :
            if isFirst : header.append(line)
        else :
            recList.append(VcfRecord(line))



def main() :

    outfp = sys.stdout

    (options,args) = getOptions()

    header=[]
    recList=[]

    isFirst=True
    for arg in args :
        processFile(arg,isFirst,header,recList)
        isFirst-False

    recList.sort(key = lambda x: (x.chrom, x.pos))

    for line in header :
        outfp.write(line)

    for vcfrec in recList :
        outfp.write(vcfrec.line)



main()