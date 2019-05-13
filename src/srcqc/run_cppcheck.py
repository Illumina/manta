#!/usr/bin/env python
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

# if cppcheck is found, run it on all c++ and return an error for *any* warning message:


import os
import sys
import subprocess


def which(searchFile) :
    """
    search the PATH for searchFile

    result should be the similar to *nix 'which' utility
    """
    for searchPath in os.environ["PATH"].split(os.pathsep):
        test=os.path.join(searchPath,searchFile)
        if os.path.isfile(test) : return test

    return None


def getCppcheckVersion() :
    """
    return the cppcheck version number, or None if version cannot be obtained.
    """

    proc=subprocess.Popen(["cppcheck","--version"],stdout=subprocess.PIPE)
    lines = proc.stdout.readlines()
    if len(lines) == 0 : return None
    return lines[0].split()[1]


def compareVersions(version1, version2):
    """
    Compare two version strings.

    Return
      < 0 if version1 is less than version2
      > 0 if version2 is less than version1
      0 if version1 == version2
    """

    def versionToIntArray(version) :
        return [int(i) for i in version.split('.')]

    # parse version components
    version1Nums = versionToIntArray(version1)
    version2Nums = versionToIntArray(version2)

    # pad versions to be the same length:
    maxlen = max(len(version1Nums), len(version2Nums))
    while len(version1Nums) < maxlen : version1Nums.append(0)
    while len(version2Nums) < maxlen : version2Nums.append(0)

    # compare
    for v1,v2 in zip(version1Nums, version2Nums) :
        if v1 < v2 : return -1
        if v1 > v2 : return 1

    return 0


def usage() :

    scriptName=os.path.basename(__file__)

    usageStr="""

    usage: %s cxx_root_directory

    run cppcheck on project c++ source code, return error for any unsupressed cppcheck issue

""" % (scriptName)

    sys.stderr.write(usageStr)
    sys.exit(2)



def main() :

    if len(sys.argv) != 2 :
        usage()

    scriptDir=os.path.dirname(sys.argv[0])

    srcRoot=os.path.abspath(sys.argv[1])

    # Update cwd so that custom cfg files will be found:
    os.chdir(scriptDir)

    cppcheck_path = which("cppcheck")
    if cppcheck_path is None :
        sys.exit(0)

    cppcheckVersion = getCppcheckVersion()

    minCppcheckVersion = "1.69"
    isOldVersion = (compareVersions(cppcheckVersion, minCppcheckVersion) < 0)
    if isOldVersion :
        sys.exit(0)

    # need to trace real path out of any symlinks so that cppcheck can find its runtime config info:
    cppcheck_path = os.path.realpath(cppcheck_path)

    checkCmd=[cppcheck_path]
    checkCmd.append("--enable=all")
    checkCmd.append("--std=c++11")
    checkCmd.append("--force")
    checkCmd.append("--verbose")
    checkCmd.append("--quiet")
    checkCmd.append("--inline-suppr")

    # manipulate the warning messages so that they look like gcc errors -- this enables IDE parsing of error location:
    checkCmd.append("--template={file}:{line}:1: error: {severity}:{message}")


    # passedByValue has been added to allow more use of c++11 shared_ptrs:
    suppressList=["unusedFunction", "unmatchedSuppression", "missingInclude", "purgedConfiguration", "passedByValue"]

    # cppcheck gives FP unitialized member errors when using delegating ctors
    suppressList.append("uninitMemberVar")

    # new items appearing on or around 1.87:
    suppressList.append("useStlAlgorithm")
    suppressList.append("assertWithSideEffect")
    suppressList.append("redundantAssignment")
    suppressList.append("variableScope")
    suppressList.append("constArgument")

    # In cppcheck versions 1.69 and lower (TODO how low?), there is a bug parsing the use of the '%' character
    # in boost::format as a regular mod operator. For these versions, an extra suppression is required.
    #
    if True :
        maxBoostFormatBugVersion = "1.69"
        isBoostFormatBugVersion = (compareVersions(cppcheckVersion, maxBoostFormatBugVersion) <= 0)
        if isBoostFormatBugVersion :
            suppressList.append("zerodivcond")

    # Version prior to 1.70 don't support current config format
    if True :
        minConfigVersion = "1.70"
        isMinConfigVersion = (compareVersions(cppcheckVersion, minConfigVersion) >= 0)
        if isMinConfigVersion :
            checkCmd.append("--library=boost_macros")

    # cppcheck v1.72 will identify lots of FP unused private method errors
    #
    # cppcheck 1.71 and 1.73 are known to not have this issue
    if True :
        minUnusedPrivateBugVersion = "1.72"
        maxUnusedPrivateBugVersion = "1.72"
        isUnusedPrivateBugVersion = (compareVersions(cppcheckVersion, minUnusedPrivateBugVersion) >= 0) and \
                                    (compareVersions(cppcheckVersion, maxUnusedPrivateBugVersion) <= 0)
        if isUnusedPrivateBugVersion :
            suppressList.append("unusedPrivateFunction")

    # cppcheck v1.71 gives FP unused struct member errors
    if True :
        minUnusedStructMemberBugVersion = "1.71"
        maxUnusedStructMemberBugVersion = "1.71"
        isUnusedStructMemberBugVersion = (compareVersions(cppcheckVersion, minUnusedStructMemberBugVersion) >= 0) and \
                                         (compareVersions(cppcheckVersion, maxUnusedStructMemberBugVersion) <= 0)
        if isUnusedStructMemberBugVersion :
            suppressList.append("unusedStructMember")

    for stype in suppressList :
        checkCmd.append("--suppress="+stype)


    # xml output is usful for getting a warnings id field, which is what you need to suppress it:
    #checkCmd.append("--xml")

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
