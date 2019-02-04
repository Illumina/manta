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

"""
This script configures the Manta SV analysis workflow
"""

import os,sys

if sys.version_info >= (3,0):
    import platform
    raise Exception("Manta does not currently support python3 (version %s detected)" % (platform.python_version()))

if sys.version_info < (2,6):
    import platform
    raise Exception("Manta requires python2 version 2.6+ (version %s detected)" % (platform.python_version()))


scriptDir=os.path.abspath(os.path.dirname(__file__))
scriptName=os.path.basename(__file__)
workflowDir=os.path.abspath(os.path.join(scriptDir,"@THIS_RELATIVE_PYTHON_LIBDIR@"))

sys.path.append(workflowDir)

from configBuildTimeInfo import workflowVersion
from mantaOptions import MantaWorkflowOptionsBase
from configureUtil import BamSetChecker, groomBamList, OptParseException, validateFixExistingFileArg
from makeRunScript import makeRunScript
from mantaWorkflow import MantaWorkflow
from workflowUtil import ensureDir



class MantaWorkflowOptions(MantaWorkflowOptionsBase) :

    def workflowDescription(self) :
        return """Version: %s

This script configures the Manta SV analysis pipeline.
You must specify a BAM or CRAM file for at least one sample.
""" % (workflowVersion)


    def addWorkflowGroupOptions(self,group) :
        group.add_option("--bam","--normalBam", type="string",dest="normalBamList",metavar="FILE", action="append",
                         help="Normal sample BAM or CRAM file. May be specified more than once, multiple inputs will be treated as each BAM file representing a different sample. [optional] (no default)")
        group.add_option("--tumorBam","--tumourBam", type="string",dest="tumorBamList",metavar="FILE", action="append",
                          help="Tumor sample BAM or CRAM file. Only up to one tumor bam file accepted. [optional] (no default)")
        group.add_option("--exome", dest="isExome", action="store_true",
                         help="Set options for WES input: turn off depth filters")
        group.add_option("--rna", dest="isRNA", action="store_true",
                         help="Set options for RNA-Seq input. Must specify exactly one bam input file")
        group.add_option("--unstrandedRNA", dest="isUnstrandedRNA", action="store_true",
                         help="Set if RNA-Seq input is unstranded: Allows splice-junctions on either strand")

        MantaWorkflowOptionsBase.addWorkflowGroupOptions(self,group)


    def addExtendedGroupOptions(self,group) :
        group.add_option("--existingAlignStatsFile",
                         dest="existingAlignStatsFile", metavar="FILE",
                         help="Pre-calculated alignment statistics file. Skips alignment stats calculation.")
        group.add_option("--useExistingChromDepths",
                         dest="useExistingChromDepths", action="store_true",
                         help="Use pre-calculated chromosome depths.")
        group.add_option("--retainTempFiles",
                         dest="isRetainTempFiles", action="store_true",
                         help="Keep all temporary files (for workflow debugging)")
        group.add_option("--generateEvidenceBam",
                         dest="isGenerateSupportBam", action="store_true",
                         help="Generate a bam of supporting reads for all SVs")
        group.add_option("--outputContig", dest="isOutputContig", action="store_true",
                         help="Output assembled contig sequences in VCF file")

        MantaWorkflowOptionsBase.addExtendedGroupOptions(self,group)


    def getOptionDefaults(self) :

        self.configScriptDir=scriptDir
        defaults=MantaWorkflowOptionsBase.getOptionDefaults(self)
        defaults.update({
            'runDir' : 'MantaWorkflow',
            'isExome' : False,
            'isRNA' : False,
            'isOutputContig' : False,
            'isUnstrandedRNA' : False,
            'existingAlignStatsFile' : None,
            'useExistingChromDepths' : False,
            'isRetainTempFiles' : False,
            'isGenerateSupportBam' : False
                          })
        return defaults


    def validateAndSanitizeOptions(self,options) :

        MantaWorkflowOptionsBase.validateAndSanitizeOptions(self,options)

        def safeLen(x) :
            if x is None : return 0
            return len(x)

        if ((safeLen(options.normalBamList) == 0) and
                (safeLen(options.tumorBamList) == 0)) :
            raise OptParseException("No normal or tumor sample alignment files specified")

        if (safeLen(options.tumorBamList) > 1) :
            raise OptParseException("Can't accept more then one tumor sample")

        if ((safeLen(options.tumorBamList) > 0) and (safeLen(options.normalBamList) > 1)) :
            raise OptParseException("Can't accept multiple normal samples for tumor subtraction")

        if options.isRNA :
            if ((safeLen(options.normalBamList) != 1) or
                    (safeLen(options.tumorBamList) != 0)) :
                raise OptParseException("RNA mode currently requires exactly one normal sample")
        else :
            if options.isUnstrandedRNA :
                raise OptParseException("Unstranded only applied for RNA inputs")

        if options.existingAlignStatsFile is not None :
            options.existingAlignStatsFile=validateFixExistingFileArg(options.existingAlignStatsFile,"existing align stats")

        groomBamList(options.normalBamList,"normal sample")
        groomBamList(options.tumorBamList, "tumor sample")

        bamSetChecker = BamSetChecker()
        if safeLen(options.normalBamList) > 0 :
            bamSetChecker.appendBams(options.normalBamList,"Normal")
        if safeLen(options.tumorBamList) > 0 :
            bamSetChecker.appendBams(options.tumorBamList,"Tumor")
        bamSetChecker.check(options.htsfileBin,
                            options.referenceFasta)



def main() :

    primarySectionName="manta"
    options,iniSections=MantaWorkflowOptions().getRunOptions(primarySectionName, version=workflowVersion)

    # We don't need to instantiate the workflow object during configuration,
    # but do it here anyway to trigger additional parameter validation:
    #
    MantaWorkflow(options)

    # generate runscript:
    #
    ensureDir(options.runDir)
    workflowScriptPath = os.path.join(options.runDir, options.workflowScriptName)

    makeRunScript(workflowScriptPath,os.path.join(workflowDir,"mantaWorkflow.py"),"MantaWorkflow",primarySectionName,iniSections)

    notefp=sys.stdout
    notefp.write("""
Successfully created workflow run script.
To execute the workflow, run the following script and set appropriate options:

%s
""" % (workflowScriptPath))


if __name__ == "__main__" :
    main()
