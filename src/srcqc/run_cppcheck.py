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

# if cppcheck is found, run it on all c++ and return an error for *any* warning message:


import os
import sys



def which(searchFile) :
    """
    search the PATH for searchFile

    result should be the similar to *nix 'which' utility
    """
    for searchPath in os.environ["PATH"].split(os.pathsep):
        test=os.path.join(searchPath,searchFile)
        if os.path.isfile(test) : return test

    return None



def usage() :

    scriptName=os.path.basename(__file__)

    usageStr="""

    usage: %s cxx_root_directory

    run cppcheck on manta c++ source code, return error for any unsupressed cppcheck issue

""" % (scriptName)

    sys.stderr.write(usageStr)
    sys.exit(2)



def main() :

    import subprocess

    if len(sys.argv) != 2 :
        usage()

    srcRoot=sys.argv[1]

    cppcheck_path = which("cppcheck")
    if cppcheck_path is None : sys.exit(0)

    # need to trace real path out of any symlinks so that cppcheck can find its runtime config info:
    cppcheck_path = os.path.realpath(cppcheck_path)

    checkCmd=[cppcheck_path]
    checkCmd.append("--enable=all")
    checkCmd.append("--std=c++03")
    checkCmd.append("--force")
    checkCmd.append("--verbose")
    checkCmd.append("--quiet")

    # manipulate the warning messages so that they look like gcc errors -- this enables IDE parsing of error location:
    checkCmd.append("--template={file}:{line}:1: error: {severity}:{message}")

    suppressList=["unusedFunction", "unmatchedSuppression", "missingInclude"]
    for stype in suppressList :
        checkCmd.append("--suppress="+stype)

    # xml output is usful for getting a warnings id field, which is what you need to suppress it:
    #checkCmd.append("--xml")

    # additional suppressed checks from starka:
    #--suppress=uninitMemberVar \
    #--suppress=unsignedLessThanZero \
    #--suppress=obsoleteFunctionsasctime \

    # this is more aggressive  and includes more FPs
    #checkCmd.append("--inconclusive")

    checkCmd.append(srcRoot)

    proc = subprocess.Popen(checkCmd, stderr=subprocess.PIPE)

    includeError="Cppcheck cannot find all the include files"

    errCount=0
    for line in proc.stderr :

        # in cppcheck 1.59 missingInclude supression appears to be broken -- this is a workaround:
        if line.find(includeError) != -1 : continue

        sys.stderr.write(line)
        errCount += 1

    if errCount != 0 : sys.exit(1)

    open("cppcheck.done", 'w').close()



if __name__ == '__main__' :
    main()
