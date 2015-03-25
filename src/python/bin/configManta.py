#!/usr/bin/env python
#
# Manta
# Copyright (c) 2013-2015 Illumina, Inc.
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
workflowDir=os.path.abspath(os.path.join(scriptDir,"@THIS_RELATIVE_PYTHON_LIBDIR@"))

version="@MANTA_FULL_VERSION@"


sys.path.append(workflowDir)

from mantaOptions import MantaWorkflowOptionsBase
from configureUtil import BamSetChecker, groomBamList, OptParseException
from makeRunScript import makeRunScript
from mantaWorkflow import MantaWorkflow
from workflowUtil import ensureDir



class MantaWorkflowOptions(MantaWorkflowOptionsBase) :

    def workflowDescription(self) :
        return """Version: %s

This script configures the Manta SV analysis pipeline.
You must specify a BAM or CRAM file for at least one sample.
""" % (version)


    def addWorkflowGroupOptions(self,group) :
        group.add_option("--normalBam", type="string",dest="normalBamList",metavar="FILE", action="append",
                         help="Normal sample BAM or CRAM file. May be specified more than once, multiple inputs will be merged. [at least one required] (no default)")
        group.add_option("--tumorBam","--tumourBam", type="string",dest="tumorBamList",metavar="FILE", action="append",
                          help="Tumor sample BAM or CRAM file. May be specified more than once, multiple inputs will be merged. [optional] (no default)")
        group.add_option("--exome", dest="isExome", action="store_true",
                         help="Set options for WES input: turn off depth filters")
        group.add_option("--rna", dest="isRNA", action="store_true",
                         help="Set options for RNA-Seq input: turn off depth filters and don't treat "
                              "anomalous reads as SV evidence when the proper-pair bit is set.")

        MantaWorkflowOptionsBase.addWorkflowGroupOptions(self,group)


    def addExtendedGroupOptions(self,group) :
        group.add_option("--useExistingAlignStats",
                         dest="useExistingAlignStats", action="store_true",
                         help="Use pre-calculated alignment statistics.")
        group.add_option("--useExistingChromDepths",
                         dest="useExistingChromDepths", action="store_true",
                         help="Use pre-calculated chromosome depths.")
        group.add_option("--candidateBins",type="int",dest="nonlocalWorkBins",metavar="candidateBins",
                         help="Provide the total number of tasks which candidate generation "
                            " will be sub-divided into. (default: %default)")

        MantaWorkflowOptionsBase.addExtendedGroupOptions(self,group)


    def getOptionDefaults(self) :

        self.configScriptDir=scriptDir
        defaults=MantaWorkflowOptionsBase.getOptionDefaults(self)
        defaults.update({
            'runDir' : 'MantaWorkflow',
            'isExome' : False,
            'isRNA' : False,
            'useExistingAlignStats' : False,
            'useExistingChromDepths' : False,
            'nonlocalWorkBins' : 256
                          })
        return defaults



    def validateAndSanitizeExistingOptions(self,options) :

        groomBamList(options.normalBamList,"normal sample")
        groomBamList(options.tumorBamList, "tumor sample")

        MantaWorkflowOptionsBase.validateAndSanitizeExistingOptions(self,options)



    def validateOptionExistence(self,options) :

        if (options.normalBamList is None) or (len(options.normalBamList) == 0) :
            raise OptParseException("No normal sample BAM files specified")

        bcheck = BamSetChecker()
        bcheck.appendBams(options.normalBamList,"Normal")
        bcheck.appendBams(options.tumorBamList,"Tumor")
        bcheck.check(options.samtoolsBin,
                     options.referenceFasta)

        MantaWorkflowOptionsBase.validateOptionExistence(self,options)




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

