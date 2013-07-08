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
Manta SV discovery workflow
"""


import os.path
import shutil
import sys

# add script path to pull in utils in same directory:
scriptDir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(scriptDir))

# add pyflow path:
# TODO: get a more robust link to the pyflow dir at config time:
pyflowDir=os.path.join(scriptDir,"pyflow")
sys.path.append(os.path.abspath(pyflowDir))

from pyflow import WorkflowRunner
from workflowUtil import checkFile, ensureDir, preJoin, which, \
                         getChromIntervals, getFastaChromOrderSize

from configureUtil import getIniSections,dumpIniSections



def getVersion() :
    return "@MANTA_VERSION@"


__version__ = getVersion()




class GenomeSegment(object) :
    """
    organizes all variables which can change
    with each genomic segment.

    The genomic segment is defined by:

    1. chromosome
    2. begin position (1-indexed closed)
    3. end position (1-indexed closed)
    4. chromosome segment (ie. bin) number (0-indexed)
    """

    def __init__(self,chrom,beginPos,endPos,binId) :
        """
        arguments are the 4 genomic interval descriptors detailed in class documentation
        """
        self.chrom = chrom
        self.bamRegion = chrom + ':' + str(beginPos) + '-' + str(endPos)
        self.binId = binId
        self.binStr = str(binId).zfill(4)
        self.id = chrom + "_" + self.binStr



def getNextGenomeSegment(params) :
    """
    generator which iterates through all genomic segments and
    returns a segmentValues object for each one.
    """
    for (chrom,beginPos,endPos,binId) in getChromIntervals(params.chromOrder,params.chromSizes,params.binSize) :
        yield GenomeSegment(chrom,beginPos,endPos,binId)



def runStats(self,taskPrefix="",dependencies=None):

    statsPath=self.paths.getStatsPath()

    cmd = [ self.params.mantaStatsBin ]
    cmd.extend(["--output-file",statsPath])
    for bamPath in self.params.normalBamList :
        cmd.extend(["--align-file",bamPath])
    for bamPath in self.params.tumorBamList :
        cmd.extend(["--tumor-align-file",bamPath])

    nextStepWait = set()
    nextStepWait.add(self.addTask(preJoin(taskPrefix,"generateStats"),cmd,dependencies=dependencies))

    return nextStepWait



def runLocusGraph(self,taskPrefix="",dependencies=None):
    """
    Create the full SV locus graph
    """

    statsPath=self.paths.getStatsPath()
    graphPath=self.paths.getGraphPath()

    graphFilename=os.path.basename(graphPath)
    tmpGraphDir=os.path.join(self.params.workDir,graphFilename+".tmpdir")
    dirTask=self.addTask(preJoin(taskPrefix,"makeTmpDir"), "mkdir -p "+tmpGraphDir, dependencies=dependencies, isForceLocal=True)

    tmpGraphFiles = []
    graphTasks = set()

    for gseg in getNextGenomeSegment(self.params) :

        tmpGraphFiles.append(os.path.join(tmpGraphDir,graphFilename+"."+gseg.id+".bin"))
        graphCmd=[ self.params.mantaGraphBin ]
        graphCmd.extend(["--output-file", tmpGraphFiles[-1]])
        graphCmd.extend(["--align-stats",statsPath])
        graphCmd.extend(["--region",gseg.bamRegion])
        for bamPath in self.params.normalBamList :
            graphCmd.extend(["--align-file",bamPath])
        for bamPath in self.params.tumorBamList :
            graphCmd.extend(["--tumor-align-file",bamPath])

        graphTaskLabel=preJoin(taskPrefix,"makeLocusGraph_"+gseg.id)
        graphTasks.add(self.addTask(graphTaskLabel,graphCmd,dependencies=dirTask))

    mergeCmd= [ self.params.mantaGraphMergeBin ]
    mergeCmd.extend(["--output-file", graphPath])
    for gfile in tmpGraphFiles :
        mergeCmd.extend(["--graph-file", gfile])

    mergeTask=self.addTask(preJoin(taskPrefix,"mergeGraph"),mergeCmd,dependencies=graphTasks)

    rmGraphTmpCmd = "rm -rf tmpGraphDir"
    #rmTask=self.addTask(preJoin(taskPrefix,"rmGraphTmp"),rmGraphTmpCmd,dependencies=mergeTask)

    nextStepWait = set()
    nextStepWait.add(mergeTask)
    return nextStepWait



def runHyGen(self, taskPrefix="", dependencies=None) :
    """
    Run hypothesis generation on each SV locus
    """

    statsPath=self.paths.getStatsPath()
    graphPath=self.paths.getGraphPath()
    hygenDir=self.paths.getHyGenDir()

    dirTask=self.addTask(preJoin(taskPrefix,"makeHyGenDir"), "mkdir -p "+ hygenDir, dependencies=dependencies, isForceLocal=True)

    isSomatic = (len(self.params.normalBamList) and len(self.params.tumorBamList))

    hygenTasks=set()
    candidateVcfPaths = []
    somaticVcfPaths = []

    for binId in range(self.params.nonlocalWorkBins) :
        binStr = str(binId).zfill(4)
        candidateVcfPaths.append(self.paths.getHyGenCandidatePath(binStr))
        if isSomatic :
            somaticVcfPaths.append(self.paths.getHyGenSomaticPath(binStr))

        hygenCmd = [ self.params.mantaHyGenBin ]
        hygenCmd.extend(["--align-stats",statsPath])
        hygenCmd.extend(["--graph-file",graphPath])
        hygenCmd.extend(["--bin-index", str(binId)])
        hygenCmd.extend(["--bin-count", str(self.params.nonlocalWorkBins)])
        hygenCmd.extend(["--ref",self.params.referenceFasta])
        hygenCmd.extend(["--candidate-output-file", candidateVcfPaths[-1]])
        if isSomatic :
            hygenCmd.extend(["--somatic-output-file", somaticVcfPaths[-1]])


        for bamPath in self.params.normalBamList :
            hygenCmd.extend(["--align-file",bamPath])
        for bamPath in self.params.tumorBamList :
            hygenCmd.extend(["--tumor-align-file",bamPath])

        hygenTaskLabel=preJoin(taskPrefix,"generateCandidateSV_"+binStr)
        hygenTasks.add(self.addTask(hygenTaskLabel,hygenCmd,dependencies=dirTask))

    nextStepWait = hygenTasks


    def getVcfSortCmd(vcfPaths, outPath) :
        cmd  = "%s -E %s " % (sys.executable,self.params.mantaSortVcf)
        cmd += " ".join(vcfPaths)
        cmd += " | %s -c > %s && %s -p vcf %s" % (self.params.bgzipBin, outPath, self.params.tabixBin, outPath)
        return cmd

    # consolidate output:
    if len(candidateVcfPaths) :
        outPath = self.paths.getSortedCandidatePath()
        candSortCmd = getVcfSortCmd(candidateVcfPaths,outPath)
        candSortLabel=preJoin(taskPrefix,"sortCandidateSV")
        nextStepWait.add(self.addTask(candSortLabel,candSortCmd,dependencies=hygenTasks))

    if len(somaticVcfPaths) :
        outPath = self.paths.getSortedSomaticPath()
        candSortCmd = getVcfSortCmd(somaticVcfPaths,outPath)
        candSortLabel=preJoin(taskPrefix,"sortSomaticSV")
        nextStepWait.add(self.addTask(candSortLabel,candSortCmd,dependencies=hygenTasks))

    return nextStepWait



class PathInfo:
    """
    object to centralize shared workflow path names
    """

    def __init__(self, params) :
        self.params = params

    def getStatsPath(self) :
        return os.path.join(self.params.workDir,"alignmentStats.txt")

    def getGraphPath(self) :
        return os.path.join(self.params.workDir,"svLocusGraph.bin")

    def getHyGenDir(self) :
        return os.path.join(self.params.workDir,"svHyGen")

    def getHyGenCandidatePath(self, binStr) :
        return os.path.join(self.getHyGenDir(),"candidateSV.%s.vcf" % (binStr))

    def getSortedCandidatePath(self) :
        return os.path.join(self.params.resultsDir,"candidateSV.vcf.gz")

    def getHyGenSomaticPath(self, binStr) :
        return os.path.join(self.getHyGenDir(),"somaticSV.%s.vcf" % (binStr))

    def getSortedSomaticPath(self) :
        return os.path.join(self.params.resultsDir,"somaticSV.vcf.gz")



class MantaWorkflow(WorkflowRunner) :
    """
    Manta SV discovery workflow
    """

    def __init__(self,params,iniSections) :
        self.params=params
        self.iniSections=iniSections

        # make sure run directory is setup:
        self.params.runDir=os.path.abspath(self.params.runDir)
        ensureDir(self.params.runDir)

        # everything that's not intended to be a final result should dump directories/files in workDir
        self.params.workDir=os.path.join(self.params.runDir,"workspace")
        ensureDir(self.params.workDir)

        # all finalized pretty results get transfered to resultsDir
        self.params.resultsDir=os.path.join(self.params.runDir,"results")
        ensureDir(self.params.resultsDir)
#         self.params.statsDir=os.path.join(self.params.resultsDir,"stats")
#         ensureDir(self.params.statsDir)
#         self.params.variantsDir=os.path.join(self.params.resultsDir,"variants")
#         ensureDir(self.params.variantsDir)
#         self.params.reportsDir=os.path.join(self.params.resultsDir,"reports")
#         ensureDir(self.params.reportsDir)

        indexRefFasta=self.params.referenceFasta+".fai"

        if self.params.referenceFasta is None:
            raise Exception("No reference fasta defined.")
        else:
            checkFile(self.params.referenceFasta,"reference fasta")
            checkFile(indexRefFasta,"reference fasta index")

        self.params.normalBamList = []
        for bam in (self.params.normalBam,) :
            if bam is None : continue
            self.params.normalBamList.append(bam)

        self.params.tumorBamList = []
        for bam in (self.params.tumorBam,) :
            if bam is None : continue
            self.params.tumorBamList.append(bam)

        # read fasta index
        (self.params.chromOrder,self.params.chromSizes) = getFastaChromOrderSize(indexRefFasta)

        # sanity check some parameter typing:
        self.params.binSize = int(self.params.binSize)
        self.params.nonlocalWorkBins = int(self.params.nonlocalWorkBins)

        self.paths = PathInfo(self.params)



    def getSuccessMessage(self) :
        "Message to be included in email for successful runs"

        msg  = "Manta workflow successfully completed.\n\n"
        msg += "\tworkflow version: %s\n" % (__version__)
        return msg



    def workflow(self) :
        self.flowLog("Initiating Manta workflow version: %s" % (__version__))

        statsTasks = runStats(self)
        graphTasks = runLocusGraph(self,dependencies=statsTasks)
        hygenTasks = runHyGen(self,dependencies=graphTasks)

