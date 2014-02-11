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
util -- simple utilities shared by workflow configurations
"""


import os.path

from optparse import OptionParser



class OptParseException(Exception):
    pass



def argToBool(x) :
    """
    convert argument of unknown type to a bool:
    """
    class FalseStrings :
        val = ("", "0", "false", "f", "no", "n", "off")

    if isinstance(x, basestring) :
        return (x.lower() not in FalseStrings.val)
    return bool(x)



def pickleConfigSections(pickleConfigFile, configSections) :
    """
    write configSections object, expected to be a hash or hashes, into a pickle file
    """
    import pickle

    pickle.dump(configSections, open(pickleConfigFile, "w"))



def getConfigSections(pickleConfigFile) :
    """
    deserialize the config file and return a hash of hashes
    """

    import pickle

    if not os.path.isfile(pickleConfigFile) : return {}

    return pickle.load(open(pickleConfigFile))



def getPrimarySectionOptions(configSections,primarySection) :

    class WorkflowOptions(object) :
        pass

    options=WorkflowOptions()
    if primarySection not in configSections : return options
    for (k,v) in configSections[primarySection].items() :
        setattr(options,k,v)

    return options



def getConfigWithPrimaryOptions(pickleConfigFile,primarySection) :
    """
    Deserialize the config pickle file and return (1) a class representing the
    options of a section specified as primary (2) a hash of hashes representing
    all sections
    """
    configSections=getConfigSections(pickleConfigFile)
    options=getPrimarySectionOptions(configSections,primarySection)

    return (options,configSections)



def dumpIniSections(iniFile,iniSections) :
    """
    convert iniSections object, expected to be a hash or hashes, into an iniFile
    """
    from ConfigParser import SafeConfigParser

    config = SafeConfigParser()
    config.optionxform=str

    def clean_value(v) :
        if v is None : return ""
        else :         return str(v)

    for section in iniSections.keys():
        config.add_section(section)
        for k,v in iniSections[section].items() :
            config.set(section,k,clean_value(v))

    configfp=open(iniFile,"w")
    config.write(configfp)
    configfp.close()



def getIniSections(iniFile) :
    """
    parse the ini iniFile and return a hash of hashes
    """
    from ConfigParser import SafeConfigParser

    if not os.path.isfile(iniFile) : return {}

    config = SafeConfigParser()
    config.optionxform=str
    config.read(iniFile)

    iniSections = {}
    for section in config.sections() :
        iniSections[section] = {}
        for (k,v) in config.items(section) :
            if v == "" : v = None
            iniSections[section][k] = v

    return iniSections



class EpilogOptionParser(OptionParser) :
    """
    This extension to OptionParser fakes the epilog feature introduced
    in versions of OptionParser after python 2.4
    """
    def __init__(self, *args, **kwargs):
        self.myepilog = None
        try:
            self.myepilog = kwargs.pop('epilog')
        except KeyError:
            pass
        OptionParser.__init__(self,*args, **kwargs)

    def print_help(self,*args,**kwargs) :
        import sys,textwrap
        OptionParser.print_help(self,*args, **kwargs)
        if self.myepilog is not None :
            sys.stdout.write("\n%s\n\n" % (textwrap.fill(self.myepilog)))



def _validateFixArgHelper(x,label,checkfunc) :
    if x is not None:
        x=os.path.abspath(x)
        if not checkfunc(x) :
            raise OptParseException("Can't find %s: '%s'" % (label,x))
    return x

def validateFixExistingDirArg(argDir,label) :
    """
    convert directory arg to absolute path and check that it exists
    """
    return _validateFixArgHelper(argDir,label,os.path.isdir)

def validateFixExistingFileArg(argFile,label) :
    """
    convert file arg to absolute path and check that it exists
    """
    return _validateFixArgHelper(argFile,label,os.path.isfile)



def checkListArgRepeats(listName,itemLabel) :
    """
    screen a list argument for repeated entries
    """
    if listName is None : return
    if len(set(listName)) != len(listName) :
        raise OptParseException("Repeated %s entries" % (itemLabel))



def assertOptionExists(arg,label) :
    if arg is None:
        raise OptParseException("No %s specified" % (label))
