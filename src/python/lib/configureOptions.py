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

import os,sys

scriptDir=os.path.abspath(os.path.dirname(__file__))
scriptName=os.path.basename(__file__)

from configureUtil import EpilogOptionParser, dumpIniSections, getIniSections, OptParseException



def noArgOrError(parser,msg) :
    if len(sys.argv) <= 1 :
        parser.print_help()
        sys.exit(2)
    else :
        parser.error(msg)



class ConfigureWorkflowOptions(object) :
    """
    This class consolidates common configuration functions
    for setting up workflows, specific configurations
    overload the indicated functions to gather/validate
    specific parameter sets for their workflow.
    """

    # override methods;

    def workflowDescription(self) :
        """
        brief description of the workflow to appear in usage
        """
        return ""

    def addWorkflowGroupOptions(self,group) :
        """
        add options to OptionsGroup object which specify
        parameters which commonly change from run to run
        """
        pass

    def addExtendedGroupOptions(self,group) :
        """
        This options are expected to change less frequently and
        should live in the ini file, they will not appear if a
        default exists.
        """
        pass

    def getOptionDefaults(self) :
        """
        Provide any option defaults. This is a good place to specify
        parameters which probably don't need to be set from their
        default. This is not a great place for anything that would
        be site-specific configuration.
        """
        return {}

    def validateAndSanitizeExistingOptions(self,options) :
        """
        validate any arguments in options which are not None, and
        optionally sanitize them. This validation step is run
        *before* writing the ini file
        """
        pass

    def validateOptionExistence(self,options) :
        """
        validate that all required options actually exist

        this method is run *after* writing template ini file
        """
        pass


    # public methods:

    def getRunOptions(self, primary_section, version = None) :
        """
        primary client code interface to the finished product.
        do not override this method

        This returns a tuple of the (1) a class holding all of the
        primary run options gathered from the primary section of the ini
        file and command-line options and (2) an inifile hash-of-hashes
        reflecting all sections of the ini file.
        """

        def updateIniSections(data,newData) :
            for k in newData.keys() :
                if k not in data : data[k] = {}
                for kk in newData[k].keys() :
                    data[k][kk] = newData[k][kk]


        # first level of options are those hard coded into the python code as defaults,
        # these have the lowest precedence:
        #
        iniSections = { primary_section : self.getOptionDefaults() }

        # next is the 'global' ini file, in the same directory as the configure
        # script:
        cmdlineScriptName=os.path.basename(sys.argv[0])
        configFileName=cmdlineScriptName+".ini"

        cmdlineScriptDir=os.path.abspath(os.path.dirname(sys.argv[0]))
        globalConfigPath=os.path.join(cmdlineScriptDir,configFileName)
        updateIniSections(iniSections,getIniSections(globalConfigPath))

        # finally there is the local ini path:
        localConfigPath=os.path.join(os.path.abspath('.'),configFileName)
        updateIniSections(iniSections,getIniSections(localConfigPath))

        parser=self._getOptionParser(iniSections[primary_section],configFileName, version=version)
        (options,args) = parser.parse_args()

        if options.isAllHelp :
            parser=self._getOptionParser(iniSections[primary_section],configFileName,True, version=version)
            parser.print_help()
            sys.exit(2)

        if len(args) : # or (len(sys.argv) == 1):
            parser.print_help()
            sys.exit(2)

        try :
            # sanitize arguments before writing defaults, check for missing arguments after:
            #
            self.validateAndSanitizeExistingOptions(options)

            # write options object back into full iniSections object:
            #
            for k,v in vars(options).iteritems() :
                if k == "isAllHelp" : continue
                if k == "isWriteConfig" : continue
                iniSections[primary_section][k] = v

            # write ini file back out if required:
            if options.isWriteConfig == True :
                dumpIniSections(configFileName,iniSections)
                sys.exit(0)

            self.validateOptionExistence(options)

        except OptParseException, e :
            noArgOrError(parser,str(e))

        return options,iniSections



    # private methods:

    def _getOptionParser(self,defaults,configFileName,isAllHelp=False, version=None) :
        from optparse import OptionGroup, SUPPRESS_HELP

        description=self.workflowDescription()+"""
Configuration will produce a workflow run script which
can execute the workflow on a single node or through
sge and resume any interrupted execution.
"""

        epilog="""Default parameters are read from global and local '%s'
files if they exist. The global ini file is searched for in the
directory containing this configuration script. The local ini files is
searched for in the current working directory. Any settings in the local
file update and take precedence over global settings in case of a repeated
entry. All current default and command-line parameters may be written to the
local ini file location with the --writeConfig option.
""" % (configFileName)

        class MyOptionParser(EpilogOptionParser) :
            def format_description(self, formatter) :
                 return self.description

        parser = MyOptionParser(description=description,epilog=epilog, version=version)

        parser.set_defaults(**defaults)

        parser.add_option("--allHelp", action="store_true",dest="isAllHelp",
                          help="show all extended/hidden options")

        group = OptionGroup(parser,"Workflow options")
        self.addWorkflowGroupOptions(group)
        parser.add_option_group(group)

        class Hack(object) : isAnyHelp=False

        class MaybeHelpOptionGroup(OptionGroup) :
            """
            This extends option group to optionally hide all help
            """
            def add_option(self,*args,**kwargs) :
                if (not isAllHelp) and \
                   ('dest' in kwargs) and \
                   (kwargs['dest'] in defaults) :
                    kwargs['help'] = SUPPRESS_HELP
                else :
                    Hack.isAnyHelp=True
                OptionGroup.add_option(self,*args, **kwargs)


        secgroup = MaybeHelpOptionGroup(parser,"Extended options",
                                        "These options are either unlikely to be reset after initial site configuration or only of interest for manta development/debugging. They will not be printed here if a default exists unless --allHelp is specified")


        def hideGroup(group) :
            group.title=group.title+" (hidden)"
            group.description=None

        self.addExtendedGroupOptions(secgroup)
        if not Hack.isAnyHelp: hideGroup(secgroup)
        parser.add_option_group(secgroup)

        def maybeHelp(key,msg) :
            if isAllHelp : return msg
            return SUPPRESS_HELP

        configgroup = OptionGroup(parser,"Config options")
        configgroup.add_option("--writeConfig", action="store_true",dest="isWriteConfig",
                               help=maybeHelp("writeConfig","Write new default configuration file based on current defaults and agruments. Defaults written to: '%s'" % (configFileName)))
        if not isAllHelp : hideGroup(configgroup)
        parser.add_option_group(configgroup)

        return parser

