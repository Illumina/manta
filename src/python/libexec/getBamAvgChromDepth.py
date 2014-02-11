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
Estimate average chromosome depth from a WGS BAM file. This will
not work correctly for exome or other targeted data.
"""

__author__ = "Chris Saunders"


import os,sys

scriptDir=os.path.abspath(os.path.dirname(__file__))
libexecDir=os.path.abspath(os.path.join(scriptDir,"@MANTA_RELATIVE_LIBEXECDIR@"))
pythonLibDir=os.path.abspath(os.path.join(scriptDir,"@MANTA_RELATIVE_PYTHON_LIBDIR@"))
sys.path.append(pythonLibDir)

from workflowUtil import checkDir, checkFile



def log(msg) :
    sys.stderr.write("INFO: " + msg + "\n")



def getOptions() :

    from optparse import OptionParser

    usage = "usage: %prog [options] > depth.txt"
    parser = OptionParser(usage=usage)

    parser.add_option("--bam", type="string",dest="bamFile",
                      help="specify bam file for depth estimation (required)")

    (options,args) = parser.parse_args()

    if len(args) != 0 :
        parser.print_help()
        sys.exit(2)

    # validate input:
    if options.bamFile is None :
        parser.print_help()
        sys.exit(2)

    checkFile(options.bamFile,"input bam")

    return (options,args)


def main() :

    import subprocess

    (options,args) = getOptions()

    checkDir(libexecDir)

    samtoolsBin = os.path.join(libexecDir,'samtools')
    checkFile(samtoolsBin,"samtools binary")

    chromData = {}
    chromList = []

    if True :
        cmd = samtoolsBin + " idxstats '%s'" % (options.bamFile)
        proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        for line in proc.stdout :
            word = line.strip().split('\t')
            if word[0] == "*" : continue

            chromData[word[0]] = (int(word[1]),int(word[2]))
            chromList.append(word[0])

    min_count = 100000

    length = 0
    count = 0
    record_count = 0

    # In a first pass, attempt to subsample the genome. If this turns out to be a tiny sample,
    # then go back and run without subsampling the bam
    #
    for type in ('subsample','all') :
        import re, signal

        length = 0
        count = 0
        record_count = 0

        # match any cigar with a series of match sequences:
        matchRex = re.compile("([0-9]+)[M=X]")

        # use "-F 4" to filter out unmapped reads
        cmd = samtoolsBin + " view -F 4"

        # use "-s 0.1" to subsample the bam records to increaase sampled read diversity
        if type == 'subsample' :
            cmd += " -s 0.1"

        cmd += " '%s'" % (options.bamFile)

        proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        for line in proc.stdout :
            record_count += 1
            word = line.strip().split('\t',7)
            isFound = False
            for mr in re.finditer(matchRex, word[5]) :
                length += int(mr.group(1))
                isFound = True
            if not isFound : continue
            count += 1
            if count >= (min_count*2) : break

        # done with subprocess:
        os.kill(proc.pid, signal.SIGINT)

        if count > min_count : break


    if count <= min_count :
        raise Exception("Unexpected read length approximation results. Observation count: " + str(count) + " Bam record count: " + str(record_count) )


    outfp=sys.stdout

    avg_length = float(length)/float(count)

    for chrom in chromList :
        if chromData[chrom][0] < avg_length : continue
        depth = chromData[chrom][1]*avg_length / float(chromData[chrom][0])
        outfp.write("%s\t%.3f\t%s\t%.3f\n" % (chrom, depth, count, avg_length))



main()
