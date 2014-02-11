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



def getBamChromInfo(samtoolsBin,bam) :
    """
    Get chromosome information from bam header

    return a map of [chrom_name]=(chrom_size,chrom_order)
    """

    import subprocess

    cmd="%s view -H '%s'" % (samtoolsBin,bam)

    info = {}
    chromIndex=0

    proc=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    for line in proc.stdout :
        if not line.startswith("@SQ") : continue
        w = line.strip().split('\t')
        if len(w) < 3 :
            chromError("Unexpected bam header for file '%s'" % (bam))

        h = {}
        for word in w[1:] :
            vals=word.split(':')
            h[vals[0]] = vals[1]

        key = h["SN"]
        size = int(h["LN"])
        if size <= 0 :
            chromError("Unexpected chromosome size '%i' in bam header for file '%s'" % (size,bam))

        info[key] = (size,chromIndex)
        chromIndex += 1

    proc.wait()
    if proc.returncode != 0 :
        chromError("Failed to pipe command: '%s'" % (cmd))

    return info



def checkChromSet(samtoolsBin,referenceFasta,bamList,bamLabel=None,isReferenceLocked=False) :
    """
    Check that chromosomes in reference and input bams are consistent

    @param samtoolsBin - samtools binary
    @param referenceFasta - samtools indexed fasta file
    @param bamList - a container of indexed bams to check for consistency
    @param bamLabel - a container of labels for each bam (default is to label BAMs by index number)
    @param isReferenceLocked - if true, then the input BAMs must contain all of the chromosomes in the reference fasta

    This function closely follows the strelka input configuration step validator
    """

    if len(bamList) == 0 : return

    if bamLabel is None :
        bamLabel = [ "index%i" % (x) for x in range(len(bamList)) ]

    assert len(bamLabel) == len(bamList)

    refChromInfo = getFastaInfo(referenceFasta)

    # first bam is used as a reference:
    chromInfo = getBamChromInfo(samtoolsBin,bamList[0])
    chroms = sorted(chromInfo.keys(),key=lambda x:chromInfo[x][1])

    # check that first bam is compatible with reference:
    for chrom in chroms :
        isError=False
        if chrom not in refChromInfo :
            isError = True
        else :
            if refChromInfo[chrom] != chromInfo[chrom][0] :
                isError = True

        if isError :
            chromError("Reference fasta and '%s' BAM file conflict on chromosome: '%s'" % (bamLabel[0],chrom))

    # optionally check that BAM contains all chromosomes in reference:
    if isReferenceLocked :
        for refChrom in refChromInfo.keys() :
            if refChrom not in chroms :
                chromError("'%s' BAM file is missing reference fasta chromosome: '%s'" % (bamLabel[0],refChrom))

    # check that other bams are compatible with first bam:
    for index in range(1,len(bamList)) :
        compareChromInfo=getBamChromInfo(samtoolsBin,bamList[index])
        for chrom in chroms:
            isError=False
            if not chrom in compareChromInfo :
                isError=True
            else :
                (ln,order) = chromInfo[chrom]
                (tln,torder) = compareChromInfo[chrom]
                if ln != tln or order != torder : isError=True

            if isError :
                chromError("'%s' and '%s' BAM files have a conflict on chromosome: '%s'" % (bamLabel[0],bamLabel[index],chrom))

            del compareChromInfo[chrom]

        # check that no chromosomes are unique to the tumor:
        for chrom in compareChromInfo.keys() :
            chromError("'%s' and '%s' BAM files have a conflict on chromosome: '%s'" % (bamLabel[0],bamLabel[index],chrom))
