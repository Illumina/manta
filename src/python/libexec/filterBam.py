#! /usr/bin/env python2
#
# Manta - Structural Variant and Indel Caller
# Copyright (c) 2013-2017 Illumina, Inc.
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
Filter the input sam file,
Only keep evidence reads supporting SVs in the candidate vcf.
"""

import sys
import re
import gzip
from os.path import isfile
from optparse import OptionParser
from glob import glob
from subprocess import call


def getOptions():
    usage = "usage: %prog [options] samtools_bin candidate_vcf input_sam filtered_sam filtered_bam"
    parser = OptionParser(usage=usage)
    (options,args) = parser.parse_args()
    if len(args) != 5 :
        parser.print_help()
        sys.exit(2)

    return (options,args)


def collect_SVs(candidateVcf, svSet):

    fpVcf = gzip.open(candidateVcf, 'rb')
    for line in fpVcf:
        if line[0] == '#':
            continue

        tokens = line.split()
        svID = tokens[2]
        svSet.add(svID)
    fpVcf.close()


def filter_sam(svSet, inputSam, filteredSam):
    fpIn = open(inputSam, 'rb')
    fpOut = open(filteredSam, 'wb')
    for line in fpIn:
        if line[0] == '@':
            fpOut.write(line)
            continue

        isSkip = True
        tokens = line[:-1].split()
        numToken = len(tokens)
        for ix in xrange(numToken):
            token = tokens[ix]
            if token[:5] == "ZM:Z:":
                newStr = token[:5]
                svItems = token[5:].split(',')
                for sv in svItems:
                    svID = sv.split("|")[0]
                    # filter out SVs not in the candidate vcf
                    if svID in svSet:
                        isSkip = False
                        newStr += sv + ','
                tokens[ix] = newStr[:-1]

        # skip the read if none of its supported SVs
        # is included in the candidate vcf
        if not(isSkip):
            for ix in xrange(numToken-1):
                fpOut.write("%s\t" % tokens[ix])
            fpOut.write("%s\n" % tokens[numToken-1])

    fpOut.close()
    fpIn.close()



if __name__=='__main__':

    # Command-line args
    (options,args) = getOptions()
    samtoolsBin = args[0]
    candidateVcf = args[1]
    inputSam = args[2]
    filteredSam = args[3]
    filteredBam = args[4]

    svSet = set([])
    collect_SVs(candidateVcf, svSet)
    filter_sam(svSet, inputSam, filteredSam)

    # convert filtered sam to bam
    call([ samtoolsBin, "view", "-h", "-b",
           "-o", filteredBam, filteredSam ])
