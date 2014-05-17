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


def getOptions() :

    from optparse import OptionParser

    usage = "usage: %prog [options] < in.vcf > out.vcf"
    description="""
reset SOMATICSCORE filter to a new value
"""
    parser = OptionParser(usage=usage)

    parser.add_option("--minSS", dest="minSS", type="int",
                      help="minimum somatic score, no default (required)")

    (opt,args) = parser.parse_args()

    if (opt.minSS is None) or (len(args) != 0) :
        parser.print_help()
        sys.exit(2)

    # validate input:
    if opt.minSS < 1:
        raise Exception("Invalid minSize value: %i" % (opt.minSS))

    return (opt,args)


def main() :

    (opt,args) = getOptions()

    msstag = "MinSomaticScore"

    infp = sys.stdin
    outfp = sys.stdout

    for line in infp :
        if line[0] == "#" :
            outfp.write(line)
            continue

        w=line.strip().split('\t')
        assert(len(w) > VCFID.INFO)

        
        val = getKeyVal(w[VCFID.INFO],"SOMATICSCORE")
        if val is None :    
            outfp.write(line)
            continue

        isPass = (int(val) >= opt.minSS)

        filters = set(w[VCFID.FILTER].split(';'))

        if len(filters) == 1 and ("PASS" in filters or "." in filters) : filters = set() 

        if isPass :
            if msstag in filters : filters.remove(msstag)
        else      : filters.add(msstag)
        
        if len(filters) == 0 :
            w[VCFID.FILTER] = "PASS"
        else :
            w[VCFID.FILTER] = ';'.join(filters) 

        sys.stdout.write('\t'.join(w) + '\n')


main()
