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

################################################################################
##
## Configuration file for the opt subfolder (external tools)
##
## author Come Raczy
##
################################################################################

message(STATUS "Creating external tools in subdirectory opt")

include ("${THIS_GLOBALS_CMAKE}")
include ("${THIS_MACROS_CMAKE}")

# convenience macro to set in both current and parent scope:
macro(superset symbol value)
    set(${symbol} "${value}")
    set(${symbol} "${value}" PARENT_SCOPE)
endmacro()

#
# htslib
#
set(HTSLIB_PREFIX "htslib-1.0")
superset(HTSLIB_DIR "${CMAKE_CURRENT_BINARY_DIR}/${HTSLIB_PREFIX}")
superset(HTSLIB_LIBRARY "${HTSLIB_DIR}/libhts.a")

add_custom_command(
    OUTPUT ${HTSLIB_LIBRARY}
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${HTSLIB_DIR}"
    COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${HTSLIB_PREFIX}.tar.bz2"
    COMMAND $(MAKE) -C "${HTSLIB_DIR}" lib-static tabix bgzip 2>| htslib.build.log
    COMMENT "Building htslib library")
add_custom_target(STARKA_HTSLIB DEPENDS "${HTSLIB_LIBRARY}")

install(PROGRAMS "${HTSLIB_DIR}/tabix" DESTINATION "${THIS_LIBEXECDIR}")
install(PROGRAMS "${HTSLIB_DIR}/bgzip" DESTINATION "${THIS_LIBEXECDIR}")

#
# samtools
#
set(SAMTOOLS_PREFIX "samtools-1.0")
superset(SAMTOOLS_DIR "${CMAKE_CURRENT_BINARY_DIR}/${SAMTOOLS_PREFIX}")
superset(SAMTOOLS_LIBRARY "${SAMTOOLS_DIR}/libbam.a")
superset(SAMTOOLS_PROG "${SAMTOOLS_DIR}/samtools")

add_custom_command(
    OUTPUT ${SAMTOOLS_PROG}
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${SAMTOOLS_DIR}"
    COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${SAMTOOLS_PREFIX}.tar.bz2"
    COMMAND ${CMAKE_COMMAND} -E chdir "${CMAKE_CURRENT_BINARY_DIR}" ln -sf "${HTSLIB_PREFIX}" htslib
    COMMAND ${CMAKE_COMMAND} -E chdir "${SAMTOOLS_DIR}" $(MAKE) all 2>| samtools.build.log
    DEPENDS "${HTSLIB_LIBRARY}"
    COMMENT "Building samtools package")
add_custom_target(MANTA_SAMTOOLS DEPENDS "${SAMTOOLS_PROG}")

install(PROGRAMS "${SAMTOOLS_PROG}" DESTINATION "${THIS_LIBEXECDIR}")


#
# pyflow
#

#make sure unpacking happens regardless of whether ${VCF_TOOLS_DIR} exists
set(PYFLOW_PREFIX "pyflow-1.1.3")
set(PYFLOW_DIR "${CMAKE_CURRENT_BINARY_DIR}/${PYFLOW_PREFIX}")
set(PYFLOW_SCRIPT "${PYFLOW_DIR}/src/pyflow.py")

add_custom_command(
    OUTPUT ${PYFLOW_SCRIPT}
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${PYFLOW_DIR}"
    COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${PYFLOW_PREFIX}.tar.bz2"
    COMMAND ${CMAKE_COMMAND} -E remove -f "${PYFLOW_DIR}/src/__init__.py"
    COMMAND ${PYTHON_EXECUTABLE} -m compileall -q "${PYFLOW_DIR}/src"
    COMMENT "Extracting pyflow")

install_python_lib_dir("${PYFLOW_DIR}/src" "${THIS_PYTHON_LIBDIR}/pyflow")

add_custom_target(STARKA_PYFLOW DEPENDS "${PYFLOW_SCRIPT}")


# tie results back to parent:
#
add_dependencies(MANTA_OPT MANTA_SAMTOOLS MANTA_TABIX MANTA_PYFLOW)


#
# hoedown
#
set(HOEDOWN_PREFIX "hoedown-1.0.2_with_allext")
set(HOEDOWN_DIR "${CMAKE_CURRENT_BINARY_DIR}/${HOEDOWN_PREFIX}")
superset(HOEDOWN_PROG "${HOEDOWN_DIR}/hoedown_allext")

add_custom_command(
    OUTPUT ${HOEDOWN_PROG}
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${HOEDOWN_DIR}"
    COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${HOEDOWN_PREFIX}.tar.bz2"
    COMMAND $(MAKE) -C "${HOEDOWN_DIR}"
    COMMENT "Building hoedown markdown parser")
add_custom_target(MANTA_MARKDOWN DEPENDS "${HOEDOWN_PROG}")

