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

#
# \author Chris Saunders
#

# final installation path:
set(INSTALL_TO_DIR "${MANTA_PYTHON_LIBDIR}")

# clean staging area for all py and pyc files:
set(PYSTAGE_DIR "${CMAKE_CURRENT_BINARY_DIR}/pystage")

file(RELATIVE_PATH MANTA_RELATIVE_LIBEXECDIR "${INSTALL_TO_DIR}" "${MANTA_LIBEXECDIR}")

include("${MANTA_MACROS_CMAKE}")
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
    DEPENDS ${STAGED_PYTHON_FILES}
    COMMAND ${PYTHON_EXECUTABLE} -m compileall -q ${PYSTAGE_DIR}
    COMMENT "Compiling python library source")

add_custom_target(MANTA_PYTHONLIB_COMPILE ALL DEPENDS ${STAGED_PYTHONC_FILES})


#
# install
#
include("${MANTA_GLOBALS_CMAKE}") # get MANTA_LIBRARY_PERMISSIONS
set(PYGLOBS "*.py" "*.pyc")
foreach(PYGLOB ${PYGLOBS})
    # Note the trialing slash in DIRECTORY argument is required:
    install (DIRECTORY "${PYSTAGE_DIR}/" DESTINATION "${INSTALL_TO_DIR}" FILE_PERMISSIONS ${MANTA_LIBRARY_PERMISSIONS} FILES_MATCHING PATTERN ${PYGLOB})
endforeach()

