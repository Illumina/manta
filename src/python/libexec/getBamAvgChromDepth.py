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
Estimate average chromosome depth from a WGS BAM file. This will
not work correctly for exome or other targeted data.
"""

__author__ = "Chris Saunders"


import os,sys

libexecDir="@MANTA_FULL_LIBEXECDIR@"
pythonLibDir="@MANTA_FULL_PYTHON_LIBDIR@"
sys.path.append(pythonLibDir)

from workflowUtil import checkDir, checkFile



def log(msg) :
    sys.stderr.write("INFO: " + msg + "\n")



def getOptions() :

    from optparse import OptionParser

    usage = "usage: %prog -bam file > depth.txt"
    parser = OptionParser(usage=usage)

    parser.add_option("--bam", type="string",dest="bamFile",
                      help="specify bam file for depth estimation (required)")

    (options,args) = parser.parse_args()

    if len(args) != 0 :
        parser.print_help()
        sys.exit(2)

    # validate input:
    if options.bamFile is None :
        parse.print_help()
        sys.exit(2)

    checkFile(bamFile,"input bam")

    return (options,args)


def main() :

    import subprocess

    (options,args) = getOptions()

    checkDir(libexecDir)

    samtoolsBin = os.path.join(libexecDir,'samtools')
    checkFile(samtoolsBin,"samtools binary")

    chromInfo = {}
    chromList = []

    if True :
        cmd = samtoolsBin + " idxstats " + options.bamFile
        proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        for line in proc.stdout :
            word = line.strip().split('\t')
            if word[0] == "*" : continue

            chromInfo[word[0]] = (word[1],word[2])
            chromList.append(word[0])


    length = 0
    count = 0

    if True :
        import re, signal

        # match any cigar with a single match sequence:
        matchRex = re.compile("^([0-9]+)M$")

        cmd = samtoolsBin + " view -F 4 -s 0.1 " + options.bamFile
        proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        for line in proc.stdout :
            word = line.strip().split('\t',7)
            mr = matchRex.match(word[5])
            if mr is None : continue
            length += int(mr.group(1))
            count += 1
            if count >= 200000 : break

        # done with subprocess:
        os.kill(proc.pid, signal.SIGINT)


    if count <= 100000 :
        raise Exception("Unexpected read length approximation results")

    avg_length = float(length)/float(count)

    for chrom in chromList :
        if chromData[chrom][0] < avg_length : continue
        depth = chromData[chrom][1]*avg_length / float(chromData[chrom][0])
        outfp.write("%s\t%.3f\t%s\t%.3f\n" % (chrom, depth, count, avg_length))



main()
