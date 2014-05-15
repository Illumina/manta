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

"""
Workflow configuration options shared by multiple
configuration scripts.
"""

import os,sys

scriptDir=os.path.abspath(os.path.dirname(__file__))
scriptName=os.path.basename(__file__)

sys.path.append(scriptDir)

from configureOptions import ConfigureWorkflowOptions
from configureUtil import assertOptionExists



def cleanLocals(locals_dict) :
    """
    When passed a locals() dictionary, clean out all of the hidden keys and return
    """

    return dict((k,v) for (k,v) in locals_dict.items() if not k.startswith("__") and k != "self")



def joinFile(*arg) :
    filePath = os.path.join(*arg)
    assert os.path.isfile(filePath)
    return filePath



class MantaWorkflowOptionsBase(ConfigureWorkflowOptions) :

    def addWorkflowGroupOptions(self,group) :
        group.add_option("--runDir", type="string",dest="runDir",
                         help="Run script and run output will be written to this directory [required] (default: %default)")

    def getOptionDefaults(self) :
        """
        Set option defaults.

        Every local variable in this method becomes part of the default hash
        """
        libexecDir=os.path.abspath(os.path.join(scriptDir,"@MANTA_RELATIVE_LIBEXECDIR@"))
        assert os.path.isdir(libexecDir)

        bgzipBin=joinFile(libexecDir,"bgzip")
        samtoolsBin=joinFile(libexecDir,"samtools")
        tabixBin=joinFile(libexecDir,"tabix")

        mantaStatsBin=joinFile(libexecDir,"GetAlignmentStats")
        mantaGraphBin=joinFile(libexecDir,"EstimateSVLoci")
        mantaGraphMergeBin=joinFile(libexecDir,"MergeSVLoci")
        mantaGraphCheckBin=joinFile(libexecDir,"CheckSVLoci")
        mantaHyGenBin=joinFile(libexecDir,"GenerateSVCandidates")
        mantaGraphStatsBin=joinFile(libexecDir,"SummarizeSVLoci")
        mantaStatsSummaryBin=joinFile(libexecDir,"SummarizeAlignmentStats")

        mantaChromDepth=joinFile(libexecDir,"getBamAvgChromDepth.py")
        mantaSortVcf=joinFile(libexecDir,"sortVcf.py")
        mantaExtraSmallVcf=joinFile(libexecDir,"extractSmallIndelCandidates.py")

        # default memory request per process-type
        #
        # where different values are provided for SGE and local runs note:
        #  1. for SGE the memory limits must be greater than the highest memory use ever
        #      expected in a production run. The consequence of exceeding this limit is a failed
        #      run.
        #   2. for localhost the memory usage should be at least above the highest mean memory
        #       use ever expected in a production run. The consequence of exceeding the mean is
        #       a slow run due to swapping.
        #
        estimateMemMb=2*1024
        mergeMemMb=4*1024
        hyGenSGEMemMb=4*1024
        hyGenLocalMemMb=2*1024

        return cleanLocals(locals())



    def validateAndSanitizeExistingOptions(self,options) :

        options.runDir=os.path.abspath(options.runDir)


    def validateOptionExistence(self,options) :

        assertOptionExists(options.runDir,"run directory")
