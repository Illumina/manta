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



class MantaWorkflowOptionsBase(ConfigureWorkflowOptions) :

    def addWorkflowGroupOptions(self,group) :
        group.add_option("--runDir", type="string",dest="runDir",
                         help="Run script and run output will be written to this directory [required] (default: %default)")

    def getOptionDefaults(self) :
        """
        Set option defaults.

        Every local variable in this method becomes part of the default hash
        """
        libexecDir="@MANTA_FULL_LIBEXECDIR@"
        assert os.path.isdir(libexecDir)

        # a hacky way to get at this, but works...
        samtoolsBin=os.path.join(libexecDir,"samtools")
        assert os.path.isfile(samtoolsBin)

        mantaStatsBin=os.path.join(libexecDir,"GetAlignmentStats")
        assert os.path.isfile(mantaStatsBin)

        mantaGraphBin=os.path.join(libexecDir,"EstimateSVLoci")
        assert os.path.isfile(mantaGraphBin)

        return cleanLocals(locals())



    def validateAndSanitizeExistingOptions(self,options) :

        options.runDir=os.path.abspath(options.runDir)


    def validateOptionExistence(self,options) :

        assertOptionExists(options.runDir,"run directory")
