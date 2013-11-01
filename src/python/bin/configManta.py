#!/usr/bin/env python
#
# Manta
# Copyright (c) 2013 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
#

"""
This script configures the Manta SV analysis workflow
"""

import os,sys

scriptDir=os.path.abspath(os.path.dirname(__file__))
scriptName=os.path.basename(__file__)
workflowDir="@MANTA_FULL_PYTHON_LIBDIR@"

version="@MANTA_FULL_VERSION@"


# quick test to see if 'make' was run
# if not os.path.exists(os.path.join(scriptDir,'..',"make.complete")) :
#     sys.stderr.write("Workflow does not appear to be fully installed. Has 'make' been run to compile sox?\n")
#     sys.exit(1)

sys.path.append(workflowDir)

from mantaOptions import MantaWorkflowOptionsBase
from configureUtil import assertOptionExists, OptParseException, validateFixExistingDirArg, validateFixExistingFileArg
from makeRunScript import makeRunScript
from mantaWorkflow import MantaWorkflow
from workflowUtil import ensureDir, isValidSampleId
from checkChromSet import checkChromSet



class MantaWorkflowOptions(MantaWorkflowOptionsBase) :

    def workflowDescription(self) :
        return """Version: %s

This script configures the Manta SV analysis pipeline.
You must specify a BAM file for at least one sample.
""" % (version)

    validAlignerModes = ["bwa","isaac"]


    def addWorkflowGroupOptions(self,group) :
        group.add_option("--normalBam", type="string",dest="normalBam",metavar="FILE",
                         help="Normal sample BAM file. [required] (no default)")
        group.add_option("--tumorBam","--tumourBam", type="string",dest="tumorBam",metavar="FILE",
                          help="Tumor sample BAM file. [optional] (no default)")
#         group.add_option("--aligner", type="string",dest="alignerMode",metavar="ALIGNER",
#                          help="Aligner type. Accepted option are {%s} [required] (no default)" % (",".join(['%s' % (x) for x in self.validAlignerModes])))
        group.add_option("--exome", dest="isExome", action="store_true",
                         help="Set options for WES input: turn off depth filters")
        group.add_option("--rna", dest="isRNA", action="store_true",
                         help="Set options for RNA-Seq input: turn off depth filters and don't treat anomalous reads as SV evidence when the proper-pair bit is set.")
        group.add_option("--referenceFasta",type="string",dest="referenceFasta",metavar="FILE",
                         help="samtools-indexed reference fasta file [required] (default: %default)")

        MantaWorkflowOptionsBase.addWorkflowGroupOptions(self,group)


    def addExtendedGroupOptions(self,group) :
        group.add_option("--useExistingAlignStats",
                         dest="useExistingAlignStats", action="store_true",
                         help="Use pre-calculated alignment statistics.")
        group.add_option("--useExistingChromDepths",
                         dest="useExistingChromDepths", action="store_true",
                         help="Use pre-calculated chromosome depths.")

        MantaWorkflowOptionsBase.addExtendedGroupOptions(self,group)


    def getOptionDefaults(self) :
        self.configScriptDir=scriptDir
        defaults=MantaWorkflowOptionsBase.getOptionDefaults(self)
        defaults.update({
            'alignerMode' : "isaac",
            'runDir' : 'MantaWorkflow',
            'isExome' : False,
            'isRNA' : False,
            'useExistingAlignStats' : False,
            'useExistingChromDepths' : False,
            'binSize' : 25000000,
            'nonlocalWorkBins' : 128
                          })
        return defaults


    def validateAndSanitizeExistingOptions(self,options) :

        options.normalBam=validateFixExistingFileArg(options.normalBam,"normal sample BAM file")
        options.tumorBam=validateFixExistingFileArg(options.tumorBam,"tumor sample BAM file")

        # check for bam index files:
        for bam in (options.tumorBam,options.normalBam) :
            if bam is None : continue
            baiFile=bam+".bai"
            if not os.path.isfile(baiFile) :
                raise OptParseException("Can't find expected BAM index file: '%s'" % (baiFile))

        # check alignerMode:
        if options.alignerMode is not None :
            options.alignerMode = options.alignerMode.lower()
            if options.alignerMode not in self.validAlignerModes :
                raise OptParseException("Invalid aligner mode: '%s'" % options.alignerMode)

        options.referenceFasta=validateFixExistingFileArg(options.referenceFasta,"reference")


        # check for reference fasta index file:
        if options.referenceFasta is not None :
            faiFile=options.referenceFasta + ".fai"
            if not os.path.isfile(faiFile) :
                raise OptParseException("Can't find expected fasta index file: '%s'" % (faiFile))

        MantaWorkflowOptionsBase.validateAndSanitizeExistingOptions(self,options)



    def validateOptionExistence(self,options) :

        assertOptionExists(options.normalBam,"normal sample BAM file")
        assertOptionExists(options.alignerMode,"aligner mode")

        assertOptionExists(options.referenceFasta,"reference fasta file")

        MantaWorkflowOptionsBase.validateOptionExistence(self,options)

        # check that the reference and the two bams are using the same set of chromosomes:
        bamList=[]
        bamLabels=[]
        if options.normalBam is not None :
            bamList.append(options.normalBam)
            bamLabels.append("Normal")

        if options.tumorBam is not None :
            bamList.append(options.tumorBam)
            bamLabels.append("Tumor")

        checkChromSet(options.samtoolsBin,
                      options.referenceFasta,
                      bamList,
                      bamLabels,
                      isReferenceLocked=True)



def main() :

    primarySectionName="manta"
    options,iniSections=MantaWorkflowOptions().getRunOptions(primarySectionName, version=version)

    # we don't need to instantiate the workflow object during configuration,
    # but this is done here to trigger additional parameter validation:
    #
    MantaWorkflow(options,iniSections)

    # generate runscript:
    #
    ensureDir(options.runDir)
    scriptFile=os.path.join(options.runDir,"runWorkflow.py")

    makeRunScript(scriptFile,os.path.join(workflowDir,"mantaWorkflow.py"),"MantaWorkflow",primarySectionName,iniSections)

    notefp=sys.stdout
    notefp.write("""
Successfully created workflow run script.
To execute the workflow, run the following script and set appropriate options:

%s
""" % (scriptFile))


if __name__ == "__main__" :
    main()

