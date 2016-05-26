#! /usr/bin/env python
#
# Manta - Structural Variant and Indel Caller
# Copyright (c) 2013-2016 Illumina, Inc.
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
Merge bams listed in a file
"""

import sys
import re
from os.path import isfile
from optparse import OptionParser
from glob import glob
from subprocess import call
from shutil import copyfile

def getOptions():
    usage = "usage: %prog [options] samtools_bin bam_mask merged_bam merged_sam bam_list_file"
    parser = OptionParser(usage=usage)
    (options,args) = parser.parse_args()

    if len(args) != 5 :
        parser.print_help()
        sys.exit(2)

    return (options,args)



if __name__=='__main__':

    # Command-line args
    (options,args) = getOptions()
    samtoolsBin = args[0]
    bamMask = args[1]
    mergedBam = args[2]
    mergedSam = args[3]
    bamListFile = args[4]

    firstBam = ""
    fileCount = 0
    fpList = open(bamListFile, 'wb')
    for bam in glob(bamMask):
        fpList.write(bam + "\n")

        if not(firstBam):
            firstBam = bam
        fileCount += 1
    fpList.close()

    if fileCount > 1:
        call([ samtoolsBin, "merge", "-b",
               bamListFile, mergedBam ])
    elif fileCount == 1:
        copyfile(firstBam, mergedBam)

    call([ samtoolsBin, "view", "-h",
           "-o", mergedSam, mergedBam ])
