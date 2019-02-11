#! /usr/bin/env python2
#
# Manta - Structural Variant and Indel Caller
# Copyright (c) 2013-2019 Illumina, Inc.
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
Sort the given bam file only if it exists
"""

import sys
from os.path import isfile
from optparse import OptionParser
from subprocess import call

def getOptions():
    usage = "usage: %prog [options] samtools_bin original_bam sorted_bam"
    parser = OptionParser(usage=usage)
    (options,args) = parser.parse_args()

    if len(args) != 3 :
        parser.print_help()
        sys.exit(2)

    return (options,args)



if __name__=='__main__':

    # Command-line args
    (options,args) = getOptions()
    samtoolsBin = args[0]
    originalBam = args[1]
    sortedBam = args[2]

    if isfile(originalBam):
        retval = call([ samtoolsBin, "sort", originalBam, "-o", sortedBam])
        if retval != 0 :
            raise Exception("Failed to sort alignment file '%s'" % (originalBam))
