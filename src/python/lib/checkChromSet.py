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

import os,sys


"""
This module contains functions to check bam header and reference consistency
"""


def chromError(msg) :
    """
    Put this here as a placeholder for custom error handling later
    """
    sys.stderr.write("\n"+"CONFIGURATION ERROR:\n"+msg+"\n\n")
    sys.exit(1)



def getFastaInfo(fasta) :
    """
    check that fai file is properly formatted (not like the GATK bundle NCBI 37 fai files)

    returns hash of chrom length
    """

    fai=fasta+".fai"
    assert os.path.isfile(fai)

    info={}

    for i,line in enumerate(open(fai)) :
        w=line.strip().split()
        if len(w) != 5 :
            msg  = "Unexpected format for line number '%i' of fasta index file: '%s'\n" % (i,fai)
            msg += "\tRe-running fasta indexing may fix the issue. To do so, run: \"samtools faidx %s\"" % (fasta)
            chromError(msg)
        info[w[0]]=int(w[1])

    return info



def getBamChromInfo(htsfileBin,bam) :
    """
    Get chromosome information from bam/cram header

    return a map of [chrom_name]=(chrom_size,chrom_order), where chrom_order is zero-indexed
    """

    import subprocess

    cmd="\"%s\" -h \"%s\"" % (htsfileBin,bam)

    info = {}
    chromIndex=0

    proc=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    for line in proc.stdout :
        if not line.startswith("@SQ") : continue
        w = line.strip().split('\t')
        if len(w) < 3 :
            chromError("Unexpected BAM/CRAM header for file '%s'" % (bam))

        h = {}
        for word in w[1:] :
            vals=word.split(':', 1)
            h[vals[0]] = vals[1]

        key = h["SN"]
        size = int(h["LN"])
        if size <= 0 :
            chromError("Unexpected chromosome size '%i' in BAM/CRAM header for file '%s'" % (size,bam))

        info[key] = (size,chromIndex)
        chromIndex += 1

    proc.wait()
    if proc.returncode != 0 :
        chromError("Failed to pipe command: '%s'" % (cmd))

    return info



def getTabixChromSet(tabixBin, tabixFile) :
    """
    Return the set of chromosomes from any tabix-indexed file
    """
    import subprocess

    chromSet = set()
    tabixCmd = [tabixBin, "-l", tabixFile]
    proc=subprocess.Popen(tabixCmd, stdout=subprocess.PIPE)
    for line in proc.stdout :
        chrom = line.strip()
        chromSet.add(chrom)

    proc.stdout.close()
    proc.wait()

    return chromSet



def ordinalStr(n) :
    """
    Given the positive integer n, return the corresponding ordinal number string
    """
    assert(n>0)

    def getSuffix(n) :
        def getOrdKey(n) :
            if n < 14 : return n
            else      : return (n % 10)

        i = getOrdKey(n)
        if   i == 1: return 'st'
        elif i == 2: return 'nd'
        elif i == 3: return 'rd'
        else :       return 'th'

    return str(n) + getSuffix(n)



def checkChromSet(htsfileBin,referenceFasta,bamList,bamLabel=None,isReferenceLocked=False) :
    """
    Check that chromosomes in reference and input bam/cram(s) are consistent

    The current requirement is very stringent. The reference and all alignment files must have
    the same set of chromosomes, and they must all be the same size. If isReferenceLocked is false
    then the reference is allowed to contain extra chromosomes not found in the alignment files.

    Within the set of alignment files, there is an additional constraint that all chromosomes are
    in the same order.

    @param htsfileBin - htsfile binary
    @param referenceFasta - samtools indexed fasta file
    @param bamList - a container of indexed bam/cram(s) to check for consistency
    @param bamLabel - a container of labels for each bam/cram file (default is to label files by index number)
    @param isReferenceLocked - if true, then the input BAM/CRAMs must contain all of the chromosomes in the reference fasta

    This function closely follows the strelka input configuration step validator
    """

    if len(bamList) == 0 : return

    if bamLabel is None :
        bamLabel = [ "index%i" % (x) for x in range(len(bamList)) ]

    assert len(bamLabel) == len(bamList)

    refChromInfo = getFastaInfo(referenceFasta)

    # first bam is used as a reference:
    chromInfo = getBamChromInfo(htsfileBin,bamList[0])
    chroms = sorted(chromInfo.keys(),key=lambda x:chromInfo[x][1])

    # check that first bam is compatible with reference:
    for chrom in chroms :
        if chrom not in refChromInfo :
            chromError("Reference genome mismatch: Reference fasta file is missing a chromosome found in the %s BAM/CRAM file: '%s'" % (bamLabel[0], chrom))
        else :
            if refChromInfo[chrom] != chromInfo[chrom][0] :
                chromError("Reference genome mismatch: The length of chromosome '%s' is %i in the reference fasta file but %i in the %s BAM/CRAM file" % (chrom, refChromInfo[chrom], chromInfo[chrom][0], bamLabel[0]))

    # optionally check that BAM contains all chromosomes in reference:
    if isReferenceLocked :
        for refChrom in refChromInfo.keys() :
            if refChrom not in chroms :
                chromError("Reference genome mismatch: %s BAM/CRAM file is missing a chromosome found in the reference fasta file: '%s'" % (bamLabel[0], refChrom))

    # check that other bams are compatible with first bam:
    for index in range(1,len(bamList)) :
        compareChromInfo=getBamChromInfo(htsfileBin,bamList[index])
        for chrom in chroms:
            if not chrom in compareChromInfo :
                chromError("Reference genome mismatch: %s BAM/CRAM file is missing a chromosome found in the %s BAM/CRAM file: '%s'" % (bamLabel[index], bamLabel[0], chrom))
            else :
                (length, order) = chromInfo[chrom]
                (compareLength, compareOrder) = compareChromInfo[chrom]
                if length != compareLength :
                    chromError("Reference genome mismatch: The length of chromosome '%s' is %i in the %s BAM/CRAM file, but %i in the %s BAM/CRAM file" % (chrom, length, bamLabel[0], compareLength, bamLabel[index]))
                if order != compareOrder :
                    chromError("Reference genome mismatch: Chromosome '%s' is ordered %s in the %s BAM/CRAM file, but %s in the %s BAM/CRAM file" % (chrom, ordinalStr(order+1), bamLabel[0], ordinalStr(compareOrder+1), bamLabel[index]))

            del compareChromInfo[chrom]

        # If all chromosomes matched, then compareChromInfo should be empty at this point.
        #
        # Check that no chromosomes are unique to the 'compare' BAM file.
        for chrom in compareChromInfo.keys() :
            chromError("Reference genome mismatch: %s BAM/CRAM file is missing a chromosome found in the %s BAM/CRAM file: '%s'" % (bamLabel[0], bamLabel[index], chrom))
