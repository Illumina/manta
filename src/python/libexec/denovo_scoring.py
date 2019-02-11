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


import sys
from os import path
from os.path import exists, abspath, dirname, basename, splitext, join


def check_genotype(probandGT, fatherGT, motherGT):
    isConsistent = False

    fatherGTItems = fatherGT.split('/')
    motherGTItems = motherGT.split('/')
    for it1 in fatherGTItems:
        for it2 in motherGTItems:
            temp = [it1, it2]
            temp.sort()
            GT = temp[0]+'/'+temp[1]
            if GT == probandGT:
                isConsistent = True
                break

    return isConsistent


def add_dq(tokens, probandIx, dq):
    for ix in xrange(9, len(tokens)):
        if (ix == probandIx):
            tokens[ix] += ":%s" % dq
        else:
            tokens[ix] += ":."


def process_vcf(vcfFile, probandID,
                fatherID, motherID):

    vcfFile = abspath(vcfFile)
    dataDir = dirname(vcfFile)
    filePrefix = splitext(basename(vcfFile))[0]
    outFile = join(dataDir, filePrefix+".de_novo.vcf")
    statsFile = join(dataDir, filePrefix+".de_novo.stats.txt")

    fpOut = open(outFile, 'wb')
    fpStats = open(statsFile, 'wb')

    countPassed = 0
    countFiltered = 0
    consistencyDict = {}

    # parser
    isFormatAdded = False
    isIxFound = False
    colNameLine = ""
    probandIx = -1
    fatherIx = -1
    motherIx = -1

    fpVcf = open(vcfFile, 'rb')
    for line in fpVcf:
        if line[0] == '#':
            if not(isFormatAdded) and (line[:8] == "##FORMAT"):
                fpOut.write("##FORMAT=<ID=DQ,Number=1,Type=Integer,Description=\"De novo quality score\">\n")
                isFormatAdded = True
            fpOut.write(line)
            colNameLine = line
            continue
        elif not(isIxFound):
            # parse format line to get the columns of proband & parents
            tokens = colNameLine.split()
            for ix in xrange(len(tokens)):
                if tokens[ix] == probandID:
                    probandIx = ix
                elif tokens[ix] == fatherID:
                    fatherIx = ix
                elif tokens[ix] == motherID:
                    motherIx = ix

            wrongID = ""
            if probandIx == -1:
                wrongID = probandID
            if fatherIx == -1:
                wrongID += (',%s' % fatherID)
            if motherIx == -1:
                wrongID += (',%s' % motherID)

            if wrongID:
                errMsg = ('The sample ID %s does not exist in the vcf.'
                          % wrongID)
                sys.stderr.write(errMsg + '\nProgram exits.')
                sys.exit(1)

        tokens = line.split()
        format = tokens[8]

        items = format.split(':')
        GTix = -1
        for ix in xrange(len(items)):
            if items[ix] == "GT":
                GTix = ix

        items = tokens[probandIx].split(':')
        probandGT = items[GTix]

        items = tokens[fatherIx].split(':')
        fatherGT = items[GTix]

        items = tokens[motherIx].split(':')
        motherGT = items[GTix]

        # add DQ to the format string
        format += ":DQ"
        isConsistent = check_genotype(probandGT, fatherGT, motherGT)
        if not(isConsistent):
            # DQ set to 60
            add_dq(tokens, probandIx, "60")

            # stats
            filter = tokens[6]
            if filter.upper() == "PASS":
                countPassed += 1
            else:
                countFiltered += 1

            GTstring = probandGT + '-' + fatherGT + '-' + motherGT
            if not(GTstring in consistencyDict):
                consistencyDict[GTstring] = 0
            consistencyDict[GTstring] += 1
        else:
            # DQ set to 0
            add_dq(tokens, probandIx, "0")

        newLine = ""
        for i in xrange(8):
            newLine += tokens[i] + "\t"
        newLine += format
        for i in xrange(9, len(tokens)):
            newLine += "\t" + tokens[i]
        fpOut.write(newLine+"\n")

    fpVcf.close()
    fpOut.close()

    fpStats.write("# of passed SVs: %s\n" % (countPassed))
    fpStats.write("# of filtered SVs: %s\n" % (countFiltered))
    fpStats.write("probandGT-fatherGT-motherGT\tcounts\n")
    genotypes = consistencyDict.keys()
    genotypes.sort()
    for gt in genotypes:
        fpStats.write("%s\t%s\n" % (gt, consistencyDict[gt]))
    fpStats.close()


if __name__=='__main__':

    usage = "denovo_scoring.py <vcf file> <proband sample ID> <father sample ID> <mother sample ID>\n"
    if len(sys.argv) <= 4:
        sys.stderr.write(usage)
        sys.exit(1)

    vcfFile = sys.argv[1]
    probandID = sys.argv[2]
    fatherID = sys.argv[3]
    motherID = sys.argv[4]

    if not(exists(vcfFile)):
        errMsg = ('The file %s does not exist.'
                  % vcfFile)
        sys.stderr.write(errMsg + '\nProgram exits.')
        sys.exit(1)

    process_vcf(vcfFile, probandID,
                fatherID, motherID)
