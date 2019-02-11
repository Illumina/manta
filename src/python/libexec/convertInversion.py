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
import gzip
from io import BufferedReader
from subprocess import check_output
from os import path
from os.path import exists


class VcfRecord:

    def __init__(self, inline):
        tokens = inline.strip().split('\t')

        self.chr = tokens[0]
        self.pos = int(tokens[1])
        self.vid = tokens[2]
        self.ref = tokens[3]
        self.alt = tokens[4]
        self.qual = tokens[5]
        self.filter = tokens[6]
        self.info = tokens[7].split(';')
        self.others = "\t".join(tokens[8:])

        # Create a dictionary for INFO
        self.infoDict ={}
        for infoItem in self.info:
            items = infoItem.split('=')
            if len(items) == 1:
                self.infoDict[items[0]] = True
            elif len(items) > 1:
                self.infoDict[items[0]] = items[1]

        self.isINV3 = False
        self.isINV5 = False
        self.mateChr = ""
        self.matePos = -1


    def checkInversion(self):

        def getMateInfo(splitChar):
            items = self.alt.split(splitChar)
            [self.mateChr, matePos] = items[1].split(':')
            self.matePos = int(matePos)

        if self.alt.startswith('['):
            getMateInfo('[')
            if self.mateChr == self.chr:
                self.isINV5 = True
        elif self.alt.endswith(']'):
            getMateInfo(']')
            if self.mateChr == self.chr:
                self.isINV3 = True


    def makeLine(self):
        infoStr = ";".join(self.info)

        self.line = "\t".join((self.chr,
                               str(self.pos),
                               self.vid,
                               self.ref,
                               self.alt,
                               self.qual,
                               self.filter,
                               infoStr,
                               self.others
                           ))+"\n"


def scanVcf(vcfFile):

    invMateDict = {}

    if vcfFile.endswith('gz'):
        gzfp = gzip.open(vcfFile, 'rb')
        fpVcf = BufferedReader(gzfp)
    else:
        fpVcf = open(vcfFile, 'rb')

    for line in fpVcf:
        if line[0] == '#':
            continue

        vcfRec = VcfRecord(line)
        vcfRec.checkInversion()
        if vcfRec.isINV3 or vcfRec.isINV5:
            if vcfRec.vid in invMateDict:
                # update mate INFO
                invMateDict[vcfRec.vid] = vcfRec.infoDict
            else:
                mateId = vcfRec.infoDict["MATEID"]
                invMateDict[mateId] = ""

    return invMateDict


def getReference(samtools, refFasta, chrom, start, end):
    region = "%s:%d-%d" % (chrom, start, end)
    samtoolsOut = check_output([samtools, "faidx", refFasta, region])
    refSeq = ""
    for seq in samtoolsOut.split('\n'):
        if not seq.startswith(">"):
            refSeq += seq

    return refSeq.upper()


def writeLines(lines):
    for line in lines:
        sys.stdout.write(line)


def convertInversions(samtools, refFasta, vcfFile, invMateDict):
    isHeaderInfoAdded = False
    isHeaderAltAdded = False
    lineBuffer = []
    bufferedChr = ""
    bufferedPos = -1

    if vcfFile.endswith('gz'):
        gzfp = gzip.open(vcfFile, 'rb')
        fpVcf = BufferedReader(gzfp)
    else:
        fpVcf = open(vcfFile, 'rb')

    for line in fpVcf:
        if line.startswith('#'):
            if (not isHeaderInfoAdded) and line.startswith("##FORMAT="):
                sys.stdout.write("##INFO=<ID=INV3,Number=0,Type=Flag,Description=\"Inversion breakends open 3' of reported location\">\n")
                sys.stdout.write("##INFO=<ID=INV5,Number=0,Type=Flag,Description=\"Inversion breakends open 5' of reported location\">\n")
                isHeaderInfoAdded = True

            if (not isHeaderAltAdded) and line.startswith("##ALT="):
                sys.stdout.write("##ALT=<ID=INV,Description=\"Inversion\">\n")
                isHeaderAltAdded = True

            sys.stdout.write(line)
            continue

        vcfRec = VcfRecord(line)

        # skip mate record
        if vcfRec.vid in invMateDict:
            continue

        vcfRec.checkInversion()
        if vcfRec.isINV3 or vcfRec.isINV5:
            if vcfRec.isINV5:
                # adjust POS for INV5
                vcfRec.pos -= 1
                vcfRec.matePos -= 1
                vcfRec.ref = getReference(samtools, refFasta,
                                          vcfRec.chr, vcfRec.pos, vcfRec.pos)

            # update manta ID
            vidSuffix = vcfRec.vid.split("MantaBND")[1]
            idx = vidSuffix.rfind(':')
            vcfRec.vid = "MantaINV%s" % vidSuffix[:idx]

            # symbolic ALT
            vcfRec.alt = "<INV>"

            # add END
            infoEndStr = "END=%d" % vcfRec.matePos

            newInfo = [infoEndStr]
            for infoItem in vcfRec.info:
                if infoItem.startswith("SVTYPE"):
                    # change SVTYPE
                    newInfo.append("SVTYPE=INV")
                    # add SVLEN
                    infoSvLenStr = "SVLEN=%d" % (vcfRec.matePos-vcfRec.pos)
                    newInfo.append(infoSvLenStr)

                elif infoItem.startswith("CIPOS"):
                    newInfo.append(infoItem)

                    # set CIEND
                    isImprecise = "IMPRECISE" in vcfRec.infoDict
                    # for imprecise calls, set CIEND to the mate breakpoint's CIPOS
                    if isImprecise:
                        mateId = vcfRec.infoDict["MATEID"]
                        mateInfoDict = invMateDict[mateId]
                        infoCiEndStr = "CIEND=%s" % (mateInfoDict["CIPOS"])
                        newInfo.append(infoCiEndStr)
                    # for precise calls, set CIEND w.r.t HOMLEN
                    else:
                        if "HOMLEN" in vcfRec.infoDict:
                            infoCiEndStr = "CIEND=-%s,0" % vcfRec.infoDict["HOMLEN"]
                            newInfo.append(infoCiEndStr)

                elif infoItem.startswith("HOMSEQ"):
                    # update HOMSEQ for INV5
                    if vcfRec.isINV5:
                        cipos = vcfRec.infoDict["CIPOS"].split(',')
                        homSeqStart = vcfRec.pos + int(cipos[0]) + 1
                        homSeqEnd = vcfRec.pos + int(cipos[1])
                        refSeq = getReference(samtools, refFasta, vcfRec.chr,
                                              homSeqStart, homSeqEnd)
                        infoHomSeqStr = "HOMSEQ=%s" % refSeq
                        newInfo.append(infoHomSeqStr)
                    else:
                        newInfo.append(infoItem)

                # skip BND-specific tags
                elif (infoItem.startswith("MATEID") or
                      infoItem.startswith("BND_DEPTH") or
                      infoItem.startswith("MATE_BND_DEPTH")):
                    continue

                # update event ID
                elif infoItem.startswith("EVENT"):
                    eidSuffix = vcfRec.infoDict["EVENT"].split("MantaBND")[1]
                    idx = vidSuffix.rfind(':')
                    infoEventStr = "EVENT=MantaINV%s" % eidSuffix[:idx]
                    newInfo.append(infoEventStr)

                # apply all other tags
                else:
                    newInfo.append(infoItem)

            # add INV3/INV5 tag
            if vcfRec.isINV3:
                newInfo.append("INV3")
            elif vcfRec.isINV5:
                newInfo.append("INV5")

            vcfRec.info = newInfo

        vcfRec.makeLine()

        # make sure the vcf is sorted in genomic order
        if (not vcfRec.chr == bufferedChr) or (vcfRec.pos > bufferedPos):
            if lineBuffer:
                writeLines(lineBuffer)

            lineBuffer = [vcfRec.line]
            bufferedChr = vcfRec.chr
            bufferedPos = vcfRec.pos
        elif vcfRec.pos < bufferedPos:
            lineBuffer.insert(0, vcfRec.line)
        else:
            lineBuffer.append(vcfRec.line)

    if lineBuffer:
        writeLines(lineBuffer)



if __name__=='__main__':

    usage = "convertInversion.py <samtools path> <reference fasta> <vcf file>\n"
    if len(sys.argv) <= 3:
        sys.stderr.write(usage)
        sys.exit(1)

    samtools = sys.argv[1]
    refFasta = sys.argv[2]
    vcfFile = sys.argv[3]

    for inputFile in [samtools, refFasta, vcfFile]:
        if not(exists(inputFile)):
            errMsg = ('File %s does not exist.'
                      % inputFile)
            sys.stderr.write(errMsg + '\nProgram exits.')
            sys.exit(1)

    invMateDict = scanVcf(vcfFile)
    convertInversions(samtools, refFasta, vcfFile, invMateDict)
