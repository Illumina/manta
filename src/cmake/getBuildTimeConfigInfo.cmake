#
# Starka
# Copyright (c) 2009-2014 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
#

# \author Chris Saunders

#
# get various build-time configuration values -- this is information not available
# at cmake configuration time
#
# requires REDIST_DIR and CONFIG_FILE
#

# generate git describe tag

set (GETGIT_CMAKE "${REDIST_DIR}/cmake-modules-c99fd3/GetGitRevisionDescription.cmake")
include ("${GETGIT_CMAKE}")
git_describe(GIT_VERSION --match "v[0-9]*")
if (NOT GIT_VERSION)
    set(GIT_VERSION "UNKNOWN")
endif ()
set(WORKFLOW_VERSION ${GIT_VERSION})
message(STATUS "Detected workflow version: ${WORKFLOW_VERSION}")
file(WRITE ${CONFIG_FILE} "WORKFLOW_VERSION\t${WORKFLOW_VERSION}\n")

#
# get build timestamp
#
# python is a cross platform way to do this without newer cmake,
# we have compile and runtime python req anyway.
#
find_package(PythonInterp QUIET)
if (NOT PYTHONINTERP_FOUND)
    message (FATAL_ERROR "No python interpreter found")
endif()

execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import datetime;print(datetime.datetime.utcnow().isoformat())"
    OUTPUT_VARIABLE BUILD_TIME
    OUTPUT_STRIP_TRAILING_WHITESPACE)

set (BUILD_TIME "${BUILD_TIME}Z")
file(APPEND ${CONFIG_FILE} "BUILD_TIME\t${BUILD_TIME}\n")
