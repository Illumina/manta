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

################################################################################
##
## Configuration file for the opt subfolder (external tools)
##
## author Come Raczy
##
################################################################################

message(STATUS "Creating external tools in subdirectory opt")

include ("${MANTA_GLOBALS_CMAKE}")

set(MANTA_REDIST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../redist")

set(SAMTOOLS "samtools-0.1.18_no_tview")
set(SAMTOOLS_DIR "${CMAKE_CURRENT_BINARY_DIR}/${SAMTOOLS}")
set(SAMTOOLS_REDIST_DIR "${MANTA_REDIST_DIR}")

set (SAMTOOLS_LIB "${SAMTOOLS_DIR}/libbam.a")

add_custom_command(OUTPUT ${SAMTOOLS_LIB}
                   COMMAND rm -rf "${SAMTOOLS_DIR}"
                   COMMAND tar -jxf "${SAMTOOLS_REDIST_DIR}/${SAMTOOLS}.tar.bz2"
                   COMMAND $(MAKE) -C "${SAMTOOLS_DIR}" all bgzip
                   COMMENT "Building the SAM/BAM library and tools")
add_custom_target(MANTA_SAMTOOLS DEPENDS "${SAMTOOLS_LIB}")

install(PROGRAMS "${SAMTOOLS_DIR}/samtools" DESTINATION "${MANTA_LIBEXECDIR}")
install(PROGRAMS "${SAMTOOLS_DIR}/bgzip"    DESTINATION "${MANTA_LIBEXECDIR}")

#make sure unpacking happens regardless of whether ${VCF_TOOLS_DIR} exists

set(PYFLOW "pyflow-1.0.0")
set(PYFLOW_DIR "${CMAKE_CURRENT_BINARY_DIR}/${PYFLOW}")
set(PYFLOW_REDIST_DIR "${MANTA_REDIST_DIR}")

add_custom_target(MANTA_PYFLOW
                  COMMAND tar -xjf "${PYFLOW_REDIST_DIR}/${PYFLOW}.tar.bz2"
                  COMMAND touch "${PYFLOW_DIR}"
                  COMMENT "Extracting pyflow")
install (DIRECTORY ${PYFLOW_DIR}/src/ DESTINATION "${MANTA_PYTHON_LIBDIR}/pyflow" FILES_MATCHING PATTERN "*.py")

# tie results back to parent:
#
set(SAMTOOLS_DIR "${SAMTOOLS_DIR}" PARENT_SCOPE)
add_dependencies(MANTA_OPT MANTA_SAMTOOLS MANTA_PYFLOW)


