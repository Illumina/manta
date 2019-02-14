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
Manta SV discovery workflow
"""


import os.path
import sys
import shutil

# add script path to pull in utils in same directory:
scriptDir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(scriptDir))

# add pyflow path:
pyflowDir=os.path.join(scriptDir,"pyflow")
sys.path.append(os.path.abspath(pyflowDir))

from checkChromSet import getTabixChromSet
from configBuildTimeInfo import workflowVersion
from configureUtil import safeSetBool
from pyflow import WorkflowRunner
from sharedWorkflow import getMkdirCmd, getMvCmd, getRmCmd, getRmdirCmd, \
                           getDepthFromAlignments
from workflowUtil import checkFile, ensureDir, preJoin, \
                         getGenomeSegmentGroups, getNextGenomeSegment, getFastaChromOrderSize, \
                         cleanPyEnv


__version__ = workflowVersion



def getOverlapCallRegions(params, genomeRegion) :
    """
    determine a set of overlapping regions between regions arguments and the callRegions bed file

    \return a list of overlapping genomic regions for calling
    """
    import subprocess

    subregions = []

    chrom = genomeRegion["chrom"]
    genomeRegionStart = genomeRegion["start"]
    genomeRegionEnd = genomeRegion["end"]
    genomeRegionStr = ("%s:%s-%s" % (chrom, genomeRegionStart, genomeRegionEnd))

    # get overlapping subregions
    tabixCmd = [params.tabixBin, params.callRegionsBed, genomeRegionStr]
    proc = subprocess.Popen(tabixCmd, stdout=subprocess.PIPE)
    for line in proc.stdout:
        # format check of the bed file
        words = line.strip().split()
        if len(words) < 3 :
            raise Exception("Unexpected format in bed file: %s\n%s" %
                            (params.callRegionsBed, line))
        callRegionStart = int(words[1])+1
        callRegionEnd = int(words[2])
        if callRegionEnd < callRegionStart :
            raise Exception("Unexpected format in bed file: %s\n%s" %
                            (params.callRegionsBed, line))

        subregionStart = max(genomeRegionStart, callRegionStart)
        subregionEnd = min(genomeRegionEnd, callRegionEnd)
        subregions.append({"chrom": chrom, "start":subregionStart, "end":subregionEnd})

    proc.stdout.close()
    proc.wait()

    return subregions



def getCallRegions(params) :
    """
    determine
    1) a set of genomic regions for calling
    2) a set of chromosomes that are completely skipped over,
       where "skipped" means that not a single base on the chrom is requested for calling

    \return a list of genomic regions for calling
    \return a set of chromLabels which are skipped
    """
    callRegionList = []
    chromIsSkipped = set()

    # when no region selections have been made:
    if ((params.genomeRegionList is None) and
            (params.callRegionsBed is None)) :
        return (callRegionList, chromIsSkipped)

    # check chromosome coverage of "regions" arguments
    chromIsSkipped = set(params.chromOrder)
    if params.genomeRegionList is not None :
        for genomeRegion in params.genomeRegionList :
            chrom = genomeRegion["chrom"]

            if chrom not in params.chromOrder:
                raise Exception("Unexpected chromosome '%s' in the argument of target regions (--region)" %
                                (chrom))

            if chrom in chromIsSkipped :
                chromIsSkipped.remove(chrom)

    if params.callRegionsBed is None :
        return (params.genomeRegionList, chromIsSkipped)

    # check chromsome coverage based on callRegions BED file
    callChromList = []
    chromIsSkipped2 = set(params.chromOrder)

    for chrom in getTabixChromSet(params.tabixBin, params.callRegionsBed) :
        if chrom not in params.chromOrder:
            raise Exception("Unexpected chromosome '%s' in the bed file of call regions %s " %
                            (chrom, params.callRegionsBed))

        callChromList.append(chrom)
        if chrom in chromIsSkipped2 :
            chromIsSkipped2.remove(chrom)

    if params.genomeRegionList is None :
        chromIsSkipped = chromIsSkipped2

        for chrom in callChromList:
            chromRegion = {"chrom":chrom, "start":1, "end":params.chromSizes[chrom]}
            callRegions = getOverlapCallRegions(params, chromRegion)
            callRegionList.extend(callRegions)
    else:
        chromIsSkipped = chromIsSkipped | chromIsSkipped2

        for genomeRegion in params.genomeRegionList:
            chrom = genomeRegion['chrom']
            if genomeRegion["start"] is None:
                genomeRegion["start"] = 1
            if genomeRegion["end"] is None:
                genomeRegion["end"] = params.chromSizes[chrom]

            subCallRegions = getOverlapCallRegions(params, genomeRegion)
            callRegionList.extend(subCallRegions)

    return (callRegionList, chromIsSkipped)



def summarizeStats(self, taskPrefix="", dependencies=None) :
    statsPath=self.paths.getStatsPath()

    summaryTasks = set()
    # summarize stats in format that's easier for human review
    cmd = [self.params.mantaStatsSummaryBin]
    cmd.extend(["--align-stats ", statsPath])
    cmd.extend(["--output-file", self.paths.getStatsSummaryPath()])
    summarizeTask = self.addTask(preJoin(taskPrefix,"summarizeStats"),cmd,dependencies=dependencies, isForceLocal=True)
    summaryTasks.add(summarizeTask)

    return summaryTasks

def copyStats(self) :
    statsPath=self.paths.getStatsPath()
    existingStatsPath=self.params.existingAlignStatsFile

    shutil.copy(existingStatsPath, statsPath)

def runStats(self,taskPrefix="",dependencies=None) :

    statsPath=self.paths.getStatsPath()
    statsFilename=os.path.basename(statsPath)

    tmpStatsDir=statsPath+".tmpdir"

    makeTmpStatsDirCmd = getMkdirCmd() + [tmpStatsDir]
    dirTask=self.addTask(preJoin(taskPrefix,"makeTmpDir"), makeTmpStatsDirCmd, dependencies=dependencies, isForceLocal=True)

    tmpStatsFiles = []
    statsTasks = set()

    for (bamIndex,bamPath) in enumerate(self.params.normalBamList + self.params.tumorBamList) :
        indexStr = str(bamIndex).zfill(3)
        tmpStatsFiles.append(os.path.join(tmpStatsDir,statsFilename+"."+ indexStr +".xml"))

        cmd = [ self.params.mantaStatsBin ]
        cmd.extend(["--ref", self.params.referenceFasta])
        cmd.extend(["--output-file",tmpStatsFiles[-1]])
        cmd.extend(["--align-file",bamPath])
        if self.params.existingAlignStatsFile:
            cmd.extend(["--default-stats-file",self.params.existingAlignStatsFile])

        statsTasks.add(self.addTask(preJoin(taskPrefix,"generateStats_"+indexStr),cmd,dependencies=dirTask))

    cmd = [ self.params.mantaMergeStatsBin ]
    cmd.extend(["--output-file",statsPath])
    for tmpStatsFile in tmpStatsFiles :
        cmd.extend(["--align-stats-file",tmpStatsFile])

    mergeTask = self.addTask(preJoin(taskPrefix,"mergeStats"),cmd,dependencies=statsTasks,isForceLocal=True)

    nextStepWait = set()
    nextStepWait.add(mergeTask)

    if not self.params.isRetainTempFiles :
        rmStatsTmpCmd = getRmdirCmd() + [tmpStatsDir]
        rmTask=self.addTask(preJoin(taskPrefix,"removeTmpDir"),rmStatsTmpCmd,dependencies=mergeTask, isForceLocal=True)

    return nextStepWait



def mantaGetDepthFromAlignments(self,taskPrefix="getChromDepth",dependencies=None):
    bamList=[]
    if len(self.params.normalBamList) :
        bamList = self.params.normalBamList
    elif len(self.params.tumorBamList) :
        bamList = self.params.tumorBamList
    else :
        return set()

    outputPath=self.paths.getChromDepth()
    return getDepthFromAlignments(self, bamList, outputPath, taskPrefix, dependencies)



def runLocusGraph(self,taskPrefix="",dependencies=None):
    """
    Create the full SV locus graph
    """

    statsPath=self.paths.getStatsPath()
    graphPath=self.paths.getGraphPath()
    graphStatsPath=self.paths.getGraphStatsPath()

    tmpGraphDir=self.paths.getTmpGraphDir()

    makeTmpGraphDirCmd = getMkdirCmd() + [tmpGraphDir]
    dirTask = self.addTask(preJoin(taskPrefix,"makeGraphTmpDir"), makeTmpGraphDirCmd, dependencies=dependencies, isForceLocal=True)

    tmpGraphFiles = []
    graphTasks = set()

    for gsegGroup in getGenomeSegmentGroups(getNextGenomeSegment(self.params)) :
        assert(len(gsegGroup) != 0)
        gid=gsegGroup[0].id
        if len(gsegGroup) > 1 :
            gid += "_to_"+gsegGroup[-1].id
        tmpGraphFiles.append(self.paths.getTmpGraphFile(gid))
        graphCmd = [ self.params.mantaGraphBin ]
        graphCmd.extend(["--output-file", tmpGraphFiles[-1]])
        graphCmd.extend(["--align-stats",statsPath])
        for gseg in gsegGroup :
            graphCmd.extend(["--region",gseg.bamRegion])
        graphCmd.extend(["--min-candidate-sv-size", self.params.minCandidateVariantSize])
        graphCmd.extend(["--min-edge-observations", self.params.minEdgeObservations])
        graphCmd.extend(["--ref",self.params.referenceFasta])
        for bamPath in self.params.normalBamList :
            graphCmd.extend(["--align-file",bamPath])
        for bamPath in self.params.tumorBamList :
            graphCmd.extend(["--tumor-align-file",bamPath])

        if self.params.isHighDepthFilter :
            graphCmd.extend(["--chrom-depth", self.paths.getChromDepth()])

        if self.params.isIgnoreAnomProperPair :
            graphCmd.append("--ignore-anom-proper-pair")
        if self.params.isRNA :
            graphCmd.append("--rna")

        graphTask=preJoin(taskPrefix,"makeLocusGraph_"+gid)
        graphTasks.add(self.addTask(graphTask,graphCmd,dependencies=dirTask,memMb=self.params.estimateMemMb))

    if len(tmpGraphFiles) == 0 :
        raise Exception("No SV Locus graphs to create. Possible target region parse error.")

    tmpGraphFileList = self.paths.getTmpGraphFileListPath()
    tmpGraphFileListTask = preJoin(taskPrefix,"mergeLocusGraphInputList")
    self.addWorkflowTask(tmpGraphFileListTask,listFileWorkflow(tmpGraphFileList,tmpGraphFiles),dependencies=graphTasks)

    mergeCmd = [ self.params.mantaGraphMergeBin ]
    mergeCmd.extend(["--output-file", graphPath])
    mergeCmd.extend(["--graph-file-list",tmpGraphFileList])
    mergeTask = self.addTask(preJoin(taskPrefix,"mergeLocusGraph"),mergeCmd,dependencies=tmpGraphFileListTask,memMb=self.params.mergeMemMb)

    # Run a separate process to rigorously check that the final graph is valid, the sv candidate generators will check as well, but
    # this makes the check much more clear:

    checkCmd = [ self.params.mantaGraphCheckBin ]
    checkCmd.extend(["--graph-file", graphPath])
    checkTask = self.addTask(preJoin(taskPrefix,"checkLocusGraph"),checkCmd,dependencies=mergeTask,memMb=self.params.mergeMemMb)

    if not self.params.isRetainTempFiles :
        rmGraphTmpCmd = getRmdirCmd() + [tmpGraphDir]
        rmTask=self.addTask(preJoin(taskPrefix,"removeTmpDir"),rmGraphTmpCmd,dependencies=mergeTask)

    graphStatsCmd  = [self.params.mantaGraphStatsBin,"--global"]
    graphStatsCmd.extend(["--graph-file",graphPath])
    graphStatsCmd.extend(["--output-file",graphStatsPath])

    graphStatsTask = self.addTask(preJoin(taskPrefix,"locusGraphStats"),graphStatsCmd,dependencies=mergeTask,memMb=self.params.mergeMemMb)

    nextStepWait = set()
    nextStepWait.add(checkTask)
    return nextStepWait



class listFileWorkflow(WorkflowRunner) :
    """
    creates a file which enumerates the values in a list, one line per item
    """

    def __init__(self2, listFile, listItems) :
        self2.listFile = listFile
        self2.listItems = listItems

    def workflow(self2) :
        fp = open(self2.listFile, "w")
        for listItem in self2.listItems :
            fp.write(listItem+"\n")



def sortEvidenceBams(self, sortBamTasks, taskPrefix="", binStr="", dependencies=None):

    for bamIdx, _ in enumerate(self.params.normalBamList + self.params.tumorBamList):
        supportBam = self.paths.getSupportBamPath(bamIdx, binStr)
        sortedBam = self.paths.getSortedSupportBamPath(bamIdx, binStr)


        # first check the existence of the supporting bam
        # then sort the bam only if it exists
        sortBamCmd = [ sys.executable, self.params.mantaSortBam,
                      self.params.samtoolsBin, supportBam, sortedBam ]

        sortBamTask = preJoin(taskPrefix, "sortEvidenceBam_%s_%s" % (binStr, bamIdx))
        sortBamTasks.add(self.addTask(sortBamTask, sortBamCmd, dependencies=dependencies))


def sortAllVcfs(self, taskPrefix="", dependencies=None) :
    """sort/prep final vcf outputs"""

    nextStepWait = set()

    def getVcfSortCmd(vcfListFile, outPath, isDiploid, isCandidate) :
        cmd  = "\"%s\" \"%s\" " % (sys.executable, self.params.mantaSortVcf)
        cmd += "-f \"%s\"" % (vcfListFile)

        # Boolean variable isCandidate is set "True" for candidateSV.vcf
        # If it is True, commandline argument "-a" is passed on to sortBam.py to print out all vcf records
        if isCandidate:
            cmd += " -a"

        # apply the ploidy filter to diploid variants
        if isDiploid:
            tempVcf = self.paths.getTempDiploidPath()
            cmd += " > \"%s\"" % (tempVcf)
            cmd += " && \"%s\" \"%s\" \"%s\"" % (sys.executable, self.params.mantaPloidyFilter, tempVcf)

        cmd += " | \"%s\" -c > \"%s\"" % (self.params.bgzipBin, outPath)

        if isDiploid:
            cmd += " && " + " ".join(getRmCmd()) + " \"%s\"" % (self.paths.getTempDiploidPath())
        return cmd

    def getVcfTabixCmd(vcfPath) :
        return [self.params.tabixBin,"-f","-p","vcf", vcfPath]


    def sortVcfs(pathList, outPath, label, isDiploid=False, isCandidate=False) :
        if len(pathList) == 0 : return set()

        # make header modifications to first vcf in list of files to be sorted:
        headerFixTask=preJoin(taskPrefix,"fixVcfHeader_"+label)
        def getHeaderFixCmd(fileName) :
            tmpName=fileName+".reheader.tmp"
            cmd  = "\"%s\" \"%s\"" % (sys.executable, self.params.vcfCmdlineSwapper)
            cmd += ' "' + " ".join(self.params.configCommandLine) + '"'
            cmd += " < \"%s\" > \"%s\"" % (fileName,tmpName)
            cmd += " && " + " ".join(getMvCmd()) +  " \"%s\" \"%s\"" % (tmpName, fileName)
            return cmd

        self.addTask(headerFixTask, getHeaderFixCmd(pathList[0]), dependencies=dependencies, isForceLocal=True)

        vcfListFile = self.paths.getVcfListPath(label)
        inputVcfTask = self.addWorkflowTask(preJoin(taskPrefix,label+"InputList"),listFileWorkflow(vcfListFile,pathList),dependencies=headerFixTask)

        sortCmd = getVcfSortCmd(vcfListFile, outPath, isDiploid, isCandidate)
        sortTask=self.addTask(preJoin(taskPrefix,"sort_"+label),sortCmd,dependencies=inputVcfTask)

        nextStepWait.add(self.addTask(preJoin(taskPrefix,"tabix_"+label),getVcfTabixCmd(outPath),dependencies=sortTask,isForceLocal=True))
        return sortTask


    candSortTask = sortVcfs(self.candidateVcfPaths,
                            self.paths.getSortedCandidatePath(),
                            "sortCandidateSV",
                            isCandidate=True)
    sortVcfs(self.diploidVcfPaths,
             self.paths.getSortedDiploidPath(),
             "sortDiploidSV",
             isDiploid=True)
    sortVcfs(self.somaticVcfPaths,
             self.paths.getSortedSomaticPath(),
             "sortSomaticSV")
    sortVcfs(self.tumorVcfPaths,
             self.paths.getSortedTumorPath(),
             "sortTumorSV")
    sortVcfs(self.rnaVcfPaths,
             self.paths.getSortedRnaPath(),
             "sortRnaSV")

    def getExtractSmallCmd(maxSize, inPath, outPath) :
        cmd  = "\"%s\" -dc \"%s\"" % (self.params.bgzipBin, inPath)
        cmd += " | \"%s\" \"%s\" --maxSize %i" % (sys.executable, self.params.mantaExtraSmallVcf, maxSize)
        cmd += " | \"%s\" -c > \"%s\"" % (self.params.bgzipBin, outPath)
        return cmd

    def extractSmall(inPath, outPath) :
        maxSize = int(self.params.minScoredVariantSize) - 1
        if maxSize < 1 : return
        smallCmd = getExtractSmallCmd(maxSize, inPath, outPath)
        smallTask=self.addTask(preJoin(taskPrefix,"extractSmallIndels"), smallCmd, dependencies=candSortTask, isForceLocal=True)
        nextStepWait.add(self.addTask(smallTask+"_tabix", getVcfTabixCmd(outPath), dependencies=smallTask, isForceLocal=True))

    extractSmall(self.paths.getSortedCandidatePath(), self.paths.getSortedCandidateSmallIndelsPath())

    return nextStepWait



def mergeEvidenceBams(self, mergeBamTasks, taskPrefix="", dependencies=None) :

    for bamIdx, bamPath in enumerate(self.params.normalBamList + self.params.tumorBamList) :
        # merge support bams
        evidenceBamFile = self.paths.getFinalSupportBamPath(bamPath, bamIdx)
        mergeCmd = [ sys.executable, self.params.mantaMergeBam,
                     self.params.samtoolsBin,
                     self.paths.getSortedSupportBamMask(bamIdx),
                     evidenceBamFile,
                     self.paths.getSupportBamListPath(bamIdx) ]

        mergeBamTask=self.addTask(preJoin(taskPrefix,"merge_evidenceBam_%s" % (bamIdx)),
                                  mergeCmd, dependencies=dependencies)

        # index the filtered bam
        ### TODO still needs to handle the case where supportBamFile does not exist
        indexCmd = [ self.params.samtoolsBin, "index", evidenceBamFile ]
        indexBamTask = self.addTask(preJoin(taskPrefix,"index_evidenceBam_%s" % (bamIdx)),
                                    indexCmd, dependencies=mergeBamTask)
        mergeBamTasks.add(indexBamTask)



def runHyGen(self, taskPrefix="", dependencies=None) :
    """
    Run hypothesis generation on each SV locus
    """

    import copy

    statsPath=self.paths.getStatsPath()
    graphPath=self.paths.getGraphPath()
    hygenDir=self.paths.getHyGenDir()

    makeHyGenDirCmd = getMkdirCmd() + [hygenDir]
    dirTask = self.addTask(preJoin(taskPrefix,"makeHyGenDir"), makeHyGenDirCmd, dependencies=dependencies, isForceLocal=True)

    isTumorNormal = (len(self.params.normalBamList) and len(self.params.tumorBamList))
    isTumorOnly = ((not isTumorNormal) and len(self.params.tumorBamList))

    hygenTasks=set()
    if self.params.isGenerateSupportBam :
        sortEvidenceBamTasks = set()

    self.candidateVcfPaths = []
    self.diploidVcfPaths = []
    self.somaticVcfPaths = []
    self.tumorVcfPaths = []
    self.rnaVcfPaths = []

    edgeRuntimeLogPaths = []

    if True:
        binId = 0
        binStr = str(binId).zfill(4)
        self.candidateVcfPaths.append(self.paths.getHyGenCandidatePath(binStr))
        if isTumorOnly :
            self.tumorVcfPaths.append(self.paths.getHyGenTumorPath(binStr))
        elif self.params.isRNA:
            self.rnaVcfPaths.append(self.paths.getHyGenRnaPath(binStr))
        else:
            self.diploidVcfPaths.append(self.paths.getHyGenDiploidPath(binStr))
            if isTumorNormal :
                self.somaticVcfPaths.append(self.paths.getHyGenSomaticPath(binStr))

        hygenCmd = [ self.params.mantaHyGenBin ]
        hygenCmd.extend(["--threads", str(self.getNCores())])
        hygenCmd.extend(["--align-stats",statsPath])
        hygenCmd.extend(["--graph-file",graphPath])
        hygenCmd.extend(["--bin-index", str(binId)])
        hygenCmd.extend(["--bin-count", "1"])
        hygenCmd.extend(["--max-edge-count", str(self.params.graphNodeMaxEdgeCount)])
        hygenCmd.extend(["--min-candidate-sv-size", self.params.minCandidateVariantSize])
        hygenCmd.extend(["--min-candidate-spanning-count", self.params.minCandidateSpanningCount])
        hygenCmd.extend(["--min-scored-sv-size", self.params.minScoredVariantSize])
        hygenCmd.extend(["--ref",self.params.referenceFasta])
        hygenCmd.extend(["--candidate-output-file", self.candidateVcfPaths[-1]])

        # tumor-only mode
        if isTumorOnly :
            hygenCmd.extend(["--tumor-output-file", self.tumorVcfPaths[-1]])
        elif self.params.isRNA:
            hygenCmd.extend(["--rna-output-file", self.rnaVcfPaths[-1]])
        else:
            hygenCmd.extend(["--diploid-output-file", self.diploidVcfPaths[-1]])
            hygenCmd.extend(["--min-qual-score", self.params.minDiploidVariantScore])
            hygenCmd.extend(["--min-pass-qual-score", self.params.minPassDiploidVariantScore])
            hygenCmd.extend(["--min-pass-gt-score", self.params.minPassDiploidGTScore])
            # tumor/normal mode
            if isTumorNormal :
                hygenCmd.extend(["--somatic-output-file", self.somaticVcfPaths[-1]])
                hygenCmd.extend(["--min-somatic-score", self.params.minSomaticScore])
                hygenCmd.extend(["--min-pass-somatic-score", self.params.minPassSomaticScore])

        # Setup remote read retrieval for insertions:
        def isEnableRemoteReadRetrieval() :
            if isTumorOnly or isTumorNormal :
                return self.params.enableRemoteReadRetrievalForInsertionsInCancerCallingModes
            else :
                return self.params.enableRemoteReadRetrievalForInsertionsInGermlineCallingModes

        if isEnableRemoteReadRetrieval() :
            hygenCmd.append("--enable-remote-read-retrieval")

        if self.params.isHighDepthFilter :
            hygenCmd.extend(["--chrom-depth", self.paths.getChromDepth()])

        edgeRuntimeLogPaths.append(self.paths.getHyGenEdgeRuntimeLogPath(binStr))
        hygenCmd.extend(["--edge-runtime-log", edgeRuntimeLogPaths[-1]])

        hygenCmd.extend(["--edge-stats-log", self.paths.getFinalEdgeStatsPath()])
        hygenCmd.extend(["--edge-stats-report", self.paths.getFinalEdgeStatsReportPath()])

        if self.params.isGenerateSupportBam :
            hygenCmd.extend(["--evidence-bam-stub", self.paths.getSupportBamStub(binStr)])

        for bamPath in self.params.normalBamList :
            hygenCmd.extend(["--align-file", bamPath])
        for bamPath in self.params.tumorBamList :
            hygenCmd.extend(["--tumor-align-file", bamPath])

        if self.params.isIgnoreAnomProperPair :
            hygenCmd.append("--ignore-anom-proper-pair")

        if self.params.useOverlapPairEvidence:
            hygenCmd.append("--use-overlapping-pair")

        if self.params.isRNA :
            hygenCmd.append("--rna")
            if self.params.isUnstrandedRNA :
                hygenCmd.append("--unstranded")

        if self.params.isOutputContig :
            hygenCmd.append("--output-contigs")

        hygenTask = preJoin(taskPrefix,"generateCandidateSV_"+binStr)
        hygenTasks.add(self.addTask(hygenTask,hygenCmd,dependencies=dirTask, nCores=self.getNCores(), memMb=self.params.hyGenMemMb))

        if self.params.isGenerateSupportBam :
            # sort evidence bams
            #
            # TODO: if the bam is large, for efficiency, consider
            # 1) filtering the bin-specific bam first w.r.t. the final candidate vcf
            # 2) then sort the bin-specific bam and merge them
            # This would require moving the filter/sort bam jobs outside the hygen loop
            sortEvidenceBams(self, sortEvidenceBamTasks, taskPrefix=taskPrefix, binStr=binStr, dependencies=hygenTask)

    nextStepWait = copy.deepcopy(hygenTasks)

    vcfTasks = sortAllVcfs(self,taskPrefix=taskPrefix,dependencies=hygenTasks)
    nextStepWait = nextStepWait.union(vcfTasks)

    if self.params.isGenerateSupportBam :
        mergeBamTasks = set()

        # merge evidence bams
        mergeEvidenceBams(self, mergeBamTasks, taskPrefix=taskPrefix,
                          dependencies=sortEvidenceBamTasks)

        nextStepWait = nextStepWait.union(mergeBamTasks)

    #
    # sort the edge runtime logs
    #
    logListFile = self.paths.getEdgeRuntimeLogListPath()
    logListTask = preJoin(taskPrefix,"sortEdgeRuntimeLogsInputList")
    self.addWorkflowTask(logListTask,listFileWorkflow(logListFile,edgeRuntimeLogPaths),dependencies=hygenTasks)

    def getEdgeLogSortCmd(logListFile, outPath) :
        cmd  = [sys.executable, self.params.mantaSortEdgeLogs,"-f", logListFile,"-o",outPath]
        return cmd

    edgeSortCmd=getEdgeLogSortCmd(logListFile,self.paths.getSortedEdgeRuntimeLogPath())
    self.addTask(preJoin(taskPrefix,"sortEdgeRuntimeLogs"), edgeSortCmd, dependencies=logListTask, isForceLocal=True)

    if not self.params.isRetainTempFiles :
        # we could delete the temp hygenDir directory here, but it is used for debug so frequently it doesn't seem worth it at present.
        # rmDirCmd = getRmdirCmd() + [hygenDir]
        # rmDirTask=self.addTask(preJoin(taskPrefix,"removeTmpDir"),rmDirCmd,dependencies=TBD_XXX_MANY)
        pass

    return nextStepWait



class PathInfo:
    """
    object to centralize shared workflow path names
    """

    def __init__(self, params) :
        self.params = params

    def getStatsPath(self) :
        return os.path.join(self.params.workDir,"alignmentStats.xml")

    def getStatsSummaryPath(self) :
        return os.path.join(self.params.statsDir,"alignmentStatsSummary.txt")

    def getChromDepth(self) :
        return os.path.join(self.params.workDir,"chromDepth.txt")

    def getGraphPath(self) :
        return os.path.join(self.params.workDir,"svLocusGraph.bin")

    def getTmpGraphDir(self) :
        return os.path.join(self.getGraphPath()+".tmpdir")

    def getTmpGraphFile(self, gid) :
        return os.path.join(self.getTmpGraphDir(),"svLocusGraph.%s.bin" % (gid))

    def getHyGenDir(self) :
        return os.path.join(self.params.workDir,"svHyGen")

    def getHyGenCandidatePath(self, binStr) :
        return os.path.join(self.getHyGenDir(),"candidateSV.%s.vcf" % (binStr))

    def getSortedCandidatePath(self) :
        return os.path.join(self.params.variantsDir,"candidateSV.vcf.gz")

    def getSortedCandidateSmallIndelsPath(self) :
        return os.path.join(self.params.variantsDir,"candidateSmallIndels.vcf.gz")

    def getHyGenDiploidPath(self, binStr) :
        return os.path.join(self.getHyGenDir(),"diploidSV.%s.vcf" % (binStr))

    def getTempDiploidPath(self) :
        return os.path.join(self.getHyGenDir(),"diploidSV.vcf.temp")

    def getSortedDiploidPath(self) :
        return os.path.join(self.params.variantsDir,"diploidSV.vcf.gz")

    def getHyGenSomaticPath(self, binStr) :
        return os.path.join(self.getHyGenDir(),"somaticSV.%s.vcf" % (binStr))

    def getHyGenTumorPath(self, binStr) :
        return os.path.join(self.getHyGenDir(),"tumorSV.%s.vcf" % (binStr))

    def getHyGenRnaPath(self, binStr) :
        return os.path.join(self.getHyGenDir(),"rnaFusion.%s.vcf" % (binStr))

    def getSortedSomaticPath(self) :
        return os.path.join(self.params.variantsDir,"somaticSV.vcf.gz")

    def getSortedTumorPath(self) :
        return os.path.join(self.params.variantsDir,"tumorSV.vcf.gz")

    def getSortedRnaPath(self) :
        return os.path.join(self.params.variantsDir,"rnaSV.vcf.gz")

    def getHyGenEdgeRuntimeLogPath(self, binStr) :
        return os.path.join(self.getHyGenDir(),"edgeRuntimeLog.%s.txt" % (binStr))

    def getSortedEdgeRuntimeLogPath(self) :
        return os.path.join(self.params.workDir,"edgeRuntimeLog.txt")

    def getFinalEdgeStatsPath(self) :
        return os.path.join(self.params.statsDir,"svCandidateGenerationStats.xml")

    def getFinalEdgeStatsReportPath(self) :
        return os.path.join(self.params.statsDir,"svCandidateGenerationStats.tsv")

    def getGraphStatsPath(self) :
        return os.path.join(self.params.statsDir,"svLocusGraphStats.tsv")


    def getSupportBamPath(self, bamIdx, binStr):
        return os.path.join(self.getHyGenDir(),
                            "evidence_%s.bam_%s.bam" % (binStr, bamIdx))

    def getSupportBamStub(self, binStr):
        return os.path.join(self.getHyGenDir(),
                            "evidence_%s" % (binStr))

    def getSortedSupportBamPath(self, bamIdx, binStr):
        return os.path.join(self.getHyGenDir(),
                            "evidence_%s.bam_%s.sorted.bam" % (binStr, bamIdx))

    def getSortedSupportBamMask(self, bamIdx):
        return os.path.join(self.getHyGenDir(),
                            "evidence_*.bam_%s.sorted.bam" % (bamIdx))

    def getFinalSupportBamPath(self, bamPath, bamIdx):
        bamPrefix = os.path.splitext(os.path.basename(bamPath))[0]
        return os.path.join(self.params.evidenceDir,
                            "evidence_%s.%s.bam" % (bamIdx, bamPrefix))

    def getSupportBamListPath(self, bamIdx):
        return os.path.join(self.getHyGenDir(),
                            "list.evidence.bam_%s.txt" % (bamIdx))

    def getTmpGraphFileListPath(self) :
        return os.path.join(self.getTmpGraphDir(),"list.svLocusGraph.txt")

    def getVcfListPath(self, label) :
        return os.path.join(self.getHyGenDir(),"list.%s.txt" % (label))

    def getEdgeRuntimeLogListPath(self) :
        return os.path.join(self.getHyGenDir(),"list.edgeRuntimeLog.txt")



class MantaWorkflow(WorkflowRunner) :
    """
    Manta SV discovery workflow
    """

    def __init__(self,params) :

        cleanPyEnv()

        self.params=params

        # normalize boolean option input:
        safeSetBool(self.params,"enableRemoteReadRetrievalForInsertionsInGermlineCallingModes")
        safeSetBool(self.params,"enableRemoteReadRetrievalForInsertionsInCancerCallingModes")
        safeSetBool(self.params,"useOverlapPairEvidence")

        # Use RNA option for minCandidate size
        if self.params.isRNA:
            self.params.minCandidateVariantSize = self.params.rnaMinCandidateVariantSize

        # format bam lists:
        if self.params.normalBamList is None : self.params.normalBamList = []
        if self.params.tumorBamList is None : self.params.tumorBamList = []

        # make sure run directory is setup:
        self.params.runDir=os.path.abspath(self.params.runDir)
        ensureDir(self.params.runDir)

        # everything that's not intended to be a final result should dump directories/files in workDir
        self.params.workDir=os.path.join(self.params.runDir,"workspace")
        ensureDir(self.params.workDir)

        # all finalized pretty results get transfered to resultsDir
        self.params.resultsDir=os.path.join(self.params.runDir,"results")
        ensureDir(self.params.resultsDir)
        self.params.statsDir=os.path.join(self.params.resultsDir,"stats")
        ensureDir(self.params.statsDir)
        self.params.variantsDir=os.path.join(self.params.resultsDir,"variants")
        ensureDir(self.params.variantsDir)
        self.params.evidenceDir=os.path.join(self.params.resultsDir,"evidence")
        ensureDir(self.params.evidenceDir)
#         self.params.reportsDir=os.path.join(self.params.resultsDir,"reports")
#         ensureDir(self.params.reportsDir)

        indexRefFasta=self.params.referenceFasta+".fai"

        if self.params.referenceFasta is None:
            raise Exception("No reference fasta defined.")
        else:
            checkFile(self.params.referenceFasta,"reference fasta")
            checkFile(indexRefFasta,"reference fasta index")

        # read fasta index
        (self.params.chromOrder,self.params.chromSizes) = getFastaChromOrderSize(indexRefFasta)
        # determine subset of chroms where we can skip calling entirely
        (self.params.callRegionList, self.params.chromIsSkipped) = getCallRegions(self.params)

        self.paths = PathInfo(self.params)

        self.params.isHighDepthFilter = (not (self.params.isExome or self.params.isRNA))
        self.params.isIgnoreAnomProperPair = (self.params.isRNA)

        # always use overlapping pairs for RNA calling
        if (self.params.isRNA) :
            self.params.useOverlapPairEvidence = True


    def setCallMemMb(self) :
        "Setup default task memory requirements"

        if self.params.callMemMbOverride is not None :
            self.params.estimateMemMb = self.params.callMemMbOverride
            self.params.hyGenMemMb = self.params.callMemMbOverride
        else :
            if self.getRunMode() == "sge" :
                self.params.hyGenMemMb = self.params.hyGenSGEMemMb
            else :
                self.params.hyGenMemMb = self.params.hyGenLocalMemMb


    def getSuccessMessage(self) :
        "Message to be included in email for successful runs"

        msg  = "Manta workflow successfully completed.\n\n"
        msg += "\tworkflow version: %s\n" % (__version__)
        return msg



    def workflow(self) :
        self.flowLog("Initiating Manta workflow version: %s" % (__version__))
        self.setCallMemMb()

        graphTaskDependencies = set()

        statsTasks = runStats(self,taskPrefix="getAlignmentStats")
        graphTaskDependencies |= statsTasks

        summarizeStats(self, dependencies=statsTasks)

        if not ((not self.params.isHighDepthFilter) or self.params.useExistingChromDepths) :
            depthTasks = mantaGetDepthFromAlignments(self)
            graphTaskDependencies |= depthTasks

        graphTasks = runLocusGraph(self,dependencies=graphTaskDependencies)

        hygenTasks = runHyGen(self,dependencies=graphTasks)
