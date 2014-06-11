#!/usr/bin/env python
#
# Manta
# Copyright (c) 2013-2014 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
#

"""
take a subset of the manta candidate vcf which can be fed into a small variant caller

select for (1) simple insert/delete combinations and (2) length <= X
"""

import os, sys
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



class VcfRecord :
    """
    simple vcf record parser
    """

    def __init__(self, line) :
        self.line = line
        w=line.strip().split('\t')
        self.chrom=w[VCFID.CHROM]
        self.pos=int(w[VCFID.POS])
        self.ref=w[VCFID.REF]
        self.alt=w[VCFID.ALT]



def getOptions() :

    from optparse import OptionParser

    usage = "usage: %prog [options] < candidate.vcf > smallIndel.vcf"
    parser = OptionParser(usage=usage)

    parser.add_option("--maxSize", dest="maxSize", type="int",
                      help="maximum indel size, no default (required)")

    (opt,args) = parser.parse_args()

    if (opt.maxSize is None) or (len(args) != 0) :
        parser.print_help()
        sys.exit(2)

    # validate input:
    if opt.maxSize < 1:
        raise Exception("Invalid maxSize value: %i" % (opt.maxSize))

    return (opt,args)



def main() :

    infp = sys.stdin
    outfp = sys.stdout

    (options,args) = getOptions()

    for line in infp :
        if line[0] == '#' :
            outfp.write(line)
            continue

        rec = VcfRecord(line)

        # remove symbolic alleles:
        if rec.alt.find("<") != -1 : continue

        # remove translocations
        if rec.alt.find("[") != -1 : continue
        if rec.alt.find("]") != -1 : continue
        if rec.alt.find(":") != -1 : continue

        # we're assume there are no multiple alts in the candidate records
        assert( rec.alt.find(",") == -1 )

        if len(rec.ref) > (options.maxSize+1) : continue
        if len(rec.alt) > (options.maxSize+1) : continue

        outfp.write(line)



main()

