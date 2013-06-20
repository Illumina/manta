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
    """

    def __init__(self,chrom,beginPos,endPos) :
        """
        arguments are the three genomic interval descriptors detailed in class documentation
        """
        self.chrom = chrom
        self.bamRegion = chrom + ':' + str(beginPos) + '-' + str(endPos)



def getNextGenomeSegment(params) :
    """
    generator which iterates through all genomic segments and
    returns a segmentValues object for each one.
    """
    for (chrom,beginPos,endPos,_) in getChromIntervals(params.chromOrder,params.chromSizes,params.binSize) :
        yield GenomeSegment(chrom,beginPos,endPos)



def runStats(self,taskPrefix="",dependencies=None):

    statsPath=self.paths.getStatsPath()

    cmd = [ self.params.mantaStatsBin ]
    cmd.extend(["--output-file",statsPath])
    for bamPath in self.params.bamList :
        cmd.extend(["--align-file",bamPath])

    nextStepWait = set()
    nextStepWait.add(self.addTask(preJoin(taskPrefix,"generateStats"),cmd,dependencies=dependencies))

    return nextStepWait



def runLocusGraph(self,taskPrefix="",dependencies=None):

    statsPath=self.paths.getStatsPath()
    graphPath=self.paths.getGraphPath()

    graphFilename=os.path.basename(graphFile)
    tmpGraphDir=os.path.join(self.params.workDir,graphFilename+".tmpdir")
    dirTask=self.addTask(preJoin(taskPrefix,"makeTmpDir"), "mkdir -p "+tmpGraphDir, dependencies=dependencies, isForceLocal=True)

    tmpGraphFiles = []
    graphTasks= []

    for gseg in getNextGenomeSegment(self.params) :

        tmpGraphFile=os.path.join(tmpGraphDir,graphFilename+gseg.bamRegion)
        graphCmd=[ self.params.mantaGraphBin ]
        graphCmd.extend(["--output-file", tmpGraphFile])
        graphCmd.extend(["--align-stats",statsPath])
        graphCmd.region(["--region",gseg.bamRegion])
        for bamPath in self.params.bamList :
            graphCmd.extend(["--align-file",bamPath])
        graphTaskLabel=preJoin(taskPrefix,"makeGraph."+gseg.bamRegion)
        graphTasks.add(self.addTask(graphTaskLabel),graphCmd,dependencies=dirTask)

    nextStepWait = set()
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

        self.params.bamList = []
        for bam in (options.tumorBam,options.normalBam) :
            if bam is None : continue
            bamList.append(bam)

        # read fasta index
        (self.params.chromOrder,self.params.chromSizes) = getFastaChromOrderSize(indexRefFasta)

        # sanity check some parameter typing:
        self.params.binSize = int(self.params.binSize)

        self.paths = PathInfo(self.params)



    def getSuccessMessage(self) :
        "Message to be included in email for successful runs"

        msg  = "Manta workflow successfully completed.\n\n"
        msg += "\tworkflow version: %s\n" % (__version__)
        return msg



    def workflow(self) :
        self.flowLog("Initiating Manta workflow version: %s" % (__version__))

        statsTasks = runStats(self)
        graphTasks = runLocusGraph(self,dependencies=statsTask)
