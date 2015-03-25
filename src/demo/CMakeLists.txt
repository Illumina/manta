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

################################################################################
##
## file CMakeLists.txt
##
################################################################################

include ("${THIS_GLOBALS_CMAKE}")

set (THIS_STRELKADEMODIR "${THIS_DEMODIR}/manta")
install(FILES "README.txt" DESTINATION ${THIS_STRELKADEMODIR})

install(FILES "runMantaWorkflowDemo.bash"  DESTINATION ${THIS_BINDIR} PERMISSIONS ${THIS_EXECUTABLE_PERMISSIONS})

add_subdirectory (data)

