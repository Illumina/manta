#!/usr/bin/env python2
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
Given a germline VCF from Manta, add/refresh the "SampleFT" FILTER on each record which does not
include a passing FORMAT/FT in at least one sample.
"""

import os, sys
import re


class VCFID :
    CHROM = 0
    POS = 1
    REF = 3
    ALT = 4
    QUAL = 5
    FILTER = 6
    INFO = 7
    FORMAT = 8



def getOptions() :

    from optparse import OptionParser

    usage = "usage: %prog < vcf > vcf_with_updated_filters"
    parser = OptionParser(usage=usage, description=__doc__)

    (options,args) = parser.parse_args()

    if (len(args) != 0) or sys.stdin.isatty() :
        parser.print_help()
        sys.exit(2)

    return (options,args)



class Constants :
    filterLabel="SampleFT"
    filterHeaderText="##FILTER=<ID=%s,Description=\"No sample passes all the sample-level filters (at the field FORMAT/FT)\">" % (filterLabel)

    nonZeroAllele = re.compile('[1-9]')



def processVariantRecordLine(outfp, line) :
    """
    Process each VCF variant record and write results to outfp stream after potentially modifying the record's
    FILTER value
    """
    w=line.strip().split('\t')
    filters=w[VCFID.FILTER].split(';')

    assert(len(filters))
    if filters[0] == "." or filters[0] == "PASS" : filters = []

    formatTags=w[VCFID.FORMAT].split(':')
    assert(len(formatTags))
    if formatTags[0] == "." : formatTags = []

    def getFilterField() :
        if len(filters) == 0 :
            return "PASS"
        else :
            return ";".join(filters)

    def outputModifiedRecord() :
        w[VCFID.FILTER] = getFilterField()
        outfp.write("\t".join(w) + "\n")

    def addFilterAndOutput() :
        if Constants.filterLabel in filters :
            outfp.write(line)
        else :
            filters.append(Constants.filterLabel)
            outputModifiedRecord()

    def removeFilterAndOutput() :
        if Constants.filterLabel not in filters :
            outfp.write(line)
        else :
            filters.remove(Constants.filterLabel)
            outputModifiedRecord()

    try :
        ftIndex = formatTags.index("FT")
    except ValueError:
        addFilterAndOutput()
        return

    isPassed=False
    for sampleIndex in range(VCFID.FORMAT+1, len(w)) :
        sampleVals=w[sampleIndex].split(':')
        ft=sampleVals[ftIndex]

        if (ft == "PASS") :
            isPassed=True
            break

    if isPassed :
        removeFilterAndOutput()
    else :
        addFilterAndOutput()




def main() :

    (_,_) = getOptions()

    infp=sys.stdin
    outfp=sys.stdout

    isFilterDescriptionFound=False
    for line in infp :

        # Scan VCF header to determine if the NoPassedVariantGTs filter description needs to be added:
        #
        if line.startswith("##") :
            if line.startswith("##FILTER") :
                if line.find(Constants.filterLabel) != -1 :
                    isFilterDescriptionFound=True
        elif line.startswith("#") :
            if not isFilterDescriptionFound :
                outfp.write(Constants.filterHeaderText + "\n")

        if line.startswith("#") :
            outfp.write(line)
            continue

        # Handle all remaining (non-header) content:
        processVariantRecordLine(outfp, line)



main()
