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
util -- simple utilities shared by bwa/gatk workflow objects
"""

__author__ = "Chris Saunders"



import os
import re



def ensureDir(d):
    """
    make directory if it doesn't already exist, raise exception if something else is in the way:
    """
    if os.path.exists(d):
        if not os.path.isdir(d) :
            raise Exception("Can't create directory: %s" % (d))
    else :
        os.makedirs(d)



def skipJoin(sep,a,b) :
    if a == "" : return b
    elif b == "" : return a
    return a+sep+b



def preJoin(a,b) :
    return skipJoin('_',a,b)



def checkFile(filename,label="") :
    if os.path.isfile(filename) : return
    if label is None : label=""
    if label != "" : label=" "+label.strip()
    raise Exception("Can't find%s file '%s'" % (label,filename) )



def which(searchFile) :
    for searchPath in os.environ["PATH"].split(":"):
        test=os.path.join(searchPath,searchFile)
        if os.path.exists(test): return test

    return None



def isValidSampleId(sampleId) :
    return re.match("^[A-Za-z0-9_-]+$", sampleId)



def getBaiFileNames(bamFile) :
    "return (picard bai filename,samtools bai filename)"
    return (bamFile[:-(len(".bam"))]+".bai",bamFile+".bai")



def javaHeapMemReqest(self,javaMb,javaMinMb=None,overheadMb=None) :
    """
    Input is the  amount of memory requested for the java heap, output is the
    amount of java heap memory you're going to actually get, and the total process memory
    (heap+overhead), to request for the task.

    If javaMinMb is not defined, it is assumed you need to full request

    If overheadMb is not defined, it is set to the global javaTaskHeapOverheadMb value

    return (javaMb,taskMb)
    """
    if javaMinMb is None : javaMinMb=javaMb
    if overheadMb is None : overheadMb=self.params.javaTaskHeapOverheadMb

    javaMb=(self.limitMemMb(javaMb+overheadMb)-overheadMb)
    if javaMb < javaMinMb :
        raise Exception("Could not provide minimum java heap memory request for task. Minimum requested: %s Available: %s" % (str(javaMinMb),str(javaMb)))
    assert (javaMb>0)
    taskMb=(javaMb+overheadMb)
    return (javaMb,taskMb)



def getFastaChromOrderSize(fastaFile) :
    """
    assuming fasta file has an fai index,
    returns
    (chromOrder,chromSizes)
    where:
    chromOrder -- list of chromosomes in fasta order
    chromSizes -- hash of chromosome sizes
    """
    faiFile=fastaFile+".fai"
    assert os.path.isfile(faiFile)

    chromOrder=[]
    chromSizes={}
    for line in open(faiFile) :
        (chrom,size)=line.strip().split("\t",2)[:2]
        chromOrder.append(chrom)
        chromSizes[chrom]=int(size)

    return (chromOrder,chromSizes)



def getChromIntervals(chromOrder,chromSizes,segmentSize) :
    """
    generate chromosome intervals no greater than segmentSize

    chromOrder - iterable object of chromosome names
    chromSizes - a hash of chrom sizes

    return chrom,start,end,chromSegment
    where start and end are formated for use with samtools
    chromSegment is 0-indexed number of segment along each chromosome
    """
    for chrom in chromOrder :
        chromSize=chromSizes[chrom]
        chromSegments=1+((chromSize-1)/segmentSize)
        segmentBaseSize=chromSize/chromSegments
        nPlusOne=chromSize%chromSegments
        start=1
        for i in xrange(chromSegments) :
            segSize=segmentBaseSize
            if i<nPlusOne : segSize += 1
            end=min(start+(segSize-1),chromSize)
            yield (chrom,start,end,i)
            start=end+1


class PathDigger(object) :
    """
    Digs into a well-defined directory structure with prefixed
    folder names to extract all files associated with
    combinations of directory names.

    This is written primarily to go through the CASAVA 1.8 output
    structure.

    #casava 1.8 fastq example:
    fqDigger=FileDigger(['Project_','Sample_'],".fastq.gz")
    """

    def __init__(self,prefixList,targetExtension=None) :
        """
        if no target extension, then list directories at the tip of the prefix list
        """
        self.prefixList=prefixList
        self.targetExtension=targetExtension


    def getNextPath(self,basePath,depth=0,ans=tuple()) :
        """
        """
        if depth < len(self.prefixList) :
            for d in os.listdir(basePath) :
                nextDir=os.path.join(basePath,d)
                if not os.path.isdir(nextDir) : continue
                if not d.startswith(self.prefixList[depth]) : continue
                value=d[len(self.prefixList[depth]):]
                for val in self.getNextPath(nextDir,depth+1,ans+tuple([value])) :
                    yield val
        else:
            if self.targetExtension is None :
                yield ans+tuple([basePath])
            else :
                for f in os.listdir(basePath) :
                    nextPath=os.path.join(basePath,f)
                    if not os.path.isfile(nextPath) : continue
                    if not f.endswith(self.targetExtension) : continue
                    yield ans+tuple([nextPath])
