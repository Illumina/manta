#!/usr/bin/env python
#
# template for quick vcf re-filtering:
#

import sys
import re

VCF_CHROM = 0
VCF_POS = 1
VCF_REF = 3
VCF_ALT = 4
VCF_QUAL = 5
VCF_FILTER = 6
VCF_INFO = 7


def getKeyVal(string,key) :
    match=re.search("%s=([^;\t]*);?" % (key) ,string)
    if match is None : return None
    return match.group(1);


for line in sys.stdin :
    w=line.strip().split('\t')
    if len(w) > VCF_FILTER and w[VCF_FILTER] == "MinSomaticScore" :
        val = getKeyVal(w[VCF_INFO],"SOMATICSCORE")
        if val is not None :    
            if val > 25 :
                w[VCF_FILTER] = "PASS"
    sys.stdout.write('\t'.join(w) + '\n')

