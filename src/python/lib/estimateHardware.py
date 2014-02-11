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
get cpu and memory capability info from linux hosts
"""


import os



class EstException(Exception):
    pass


def getNodeRealCoreCount() :
    """
    Only works on Linux systems
    Reads /proc/cpuinfo and tries to determine the real number of cpu cores

    adapted from B Sickler
    """
    cpuinfo = '/proc/cpuinfo'
    if not os.path.isfile(cpuinfo):
        raise EstException("Can't read processor information from '%s'" % (cpuinfo))

    physical_core_ids = set()
    cpu_cores = set()

    for line in open(cpuinfo):
        l_array = [item.strip() for item in line.strip().split(': ')]
        if len(l_array) < 2: continue
        cpuid = l_array[0]
        if cpuid == 'cpu cores':
            cpu_cores.add(int(l_array[1]))
        if cpuid == 'physical id':
            physical_core_ids.add(l_array[1])

    if len(physical_core_ids) == 0 :
        raise EstException("No 'physical id' key found in '%s'" % (cpuinfo))
    if len(cpu_cores) == 0 :
        raise EstException("No 'cpu cores' key found in '%s'" % (cpuinfo))

    return ( len(physical_core_ids) * cpu_cores.pop() )




def getNodeHyperthreadCoreCount():
    """
    return the number of hyperthread cores on this host

    Taken from R Kelley's function in IsisWorkflow
    """

    #
    # We're being called from a qsub,
    # respect the core count that it was
    # trying to assign
    #
    #if "NSLOTS" in os.environ:
    #    return int(os.environ["NSLOTS"])

    #
    # User didn't specify a core limit and neither
    # did SGE. It's open season on CPU cores.
    #
    cname="/proc/cpuinfo"
    if not os.path.isfile(cname):
        raise EstException("Can't read processor information from %s" % (cname))

    cpuCount = 0
    for line in open(cname):
        if line.startswith("processor"): cpuCount += 1

    if cpuCount == 0: raise EstException("Can't estimate processor core count from %s" % (cname))

    return cpuCount



def getNodeMemMb():
    """
    return total memory in Mbytes

    Taken from R Kelley's function in IsisWorkflow
    """
    #
    # Can only get this from /proc/meminfo
    #
    mname="/proc/meminfo"
    if not os.path.isfile(mname):
        raise EstException("Can't read memory information from %s" % (mname))

    line = open(mname).readline()
    splat = line.rstrip().split()
    if len(splat) != 3:
        raise EstException("Unexpected format in %s" % (mname))

    try:
        memMb = 1+(int(splat[1])-1)/1024
    except:
        raise EstException("Unexpected format in %s" % (mname))

    return memMb


