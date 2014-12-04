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

#
# \author Chris Saunders
#

# final installation path:
set(INSTALL_TO_DIR "${THIS_PYTHON_LIBDIR}")

# clean staging area for all py and pyc files:
set(PYSTAGE_DIR "${CMAKE_CURRENT_BINARY_DIR}/pystage")

file(RELATIVE_PATH THIS_RELATIVE_LIBEXECDIR "${INSTALL_TO_DIR}" "${THIS_LIBEXECDIR}")

include("${THIS_MACROS_CMAKE}")
configure_files("${CMAKE_CURRENT_SOURCE_DIR}" "${PYSTAGE_DIR}" "*.py")


#
# compile all py to pyc (for build-time error-checking and faster run-time response)
#
file(GLOB STAGED_PYTHON_FILES "${PYSTAGE_DIR}/*.py")
foreach(PYFILE ${STAGED_PYTHON_FILES})
    list(APPEND STAGED_PYTHONC_FILES "${PYFILE}c")
endforeach()

file(GLOB STAGED_PYTHON_FILES "${PYSTAGE_DIR}/*.py")
add_custom_command(
    OUTPUT ${STAGED_PYTHONC_FILES}
    COMMAND ${PYTHON_EXECUTABLE} -m compileall -q ${PYSTAGE_DIR}
    COMMENT "Compiling python library source")

set(THIS_PYTHONLIB_COMPILE "${THIS_PROJECT_NAME}_pythonlib_compile")
add_custom_target(${THIS_PYTHONLIB_COMPILE} ALL DEPENDS ${STAGED_PYTHONC_FILES})


#
# install
#
include("${THIS_MACROS_CMAKE}")
install_python_lib_dir("${PYSTAGE_DIR}" "${INSTALL_TO_DIR}")

