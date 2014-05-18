#!/usr/bin/env python
#
# template for quick vcf re-filtering:
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
    SAMPLE = 8



def main() :

    infp = sys.stdin
    outfp = sys.stdout

    for line in infp :
        if line[0] == "#" :
            outfp.write(line)
            continue

        w=line.strip().split('\t')
        assert(len(w) > VCFID.SAMPLE+1)

        x=w[VCFID.SAMPLE+1].split(':')
        normalAltPairs=int(x[0].split(',')[1])

        val = getKeyVal(w[VCFID.INFO],"SOMATICSCORE")
        if val is not None and int(val) <= 20 and normalAltPairs != 0 : continue

        outfp.write(line)

main()
