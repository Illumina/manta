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
from workflowUtil import checkFile, ensureDir, preJoin, which

from configureUtil import getIniSections,dumpIniSections



def getVersion() :
    return "@MANTA_VERSION@"


__version__ = getVersion()



class PathInfo:
    """
    object to centralize shared workflow path names
    """

    def __init__(self, params) :
        self.params = params



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
        self.params.statsDir=os.path.join(self.params.resultsDir,"stats")
        ensureDir(self.params.statsDir)
        self.params.variantsDir=os.path.join(self.params.resultsDir,"variants")
        ensureDir(self.params.variantsDir)
        self.params.reportsDir=os.path.join(self.params.resultsDir,"reports")
        ensureDir(self.params.reportsDir)

        if self.params.referenceFasta is None:
            raise Exception("No reference fasta defined.")
        else:
            checkFile(self.params.referenceFasta,"reference fasta")
            checkFile(self.params.referenceFasta+".fai","reference fasta index")

        self.paths = PathInfo(self.params)



    def getSuccessMessage(self) :
        "Message to be included in email for successful runs"

        msg  = "Manta workflow successfully completed.\n\n"
        msg += "\tworkflow version: %s\n" % (__version__)
        return msg



    def workflow(self) :
        self.flowLog("Initiating Manta workflow version: %s" % (__version__))

