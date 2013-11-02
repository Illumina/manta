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

set(INSTALL_TO_DIR "${MANTA_PYTHON_LIBDIR}")

file(RELATIVE_PATH MANTA_RELATIVE_LIBEXECDIR "${INSTALL_TO_DIR}" "${MANTA_LIBEXECDIR}")

include("${MANTA_MACROS_CMAKE}")
configure_files("${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/pystage" "*.py")

add_custom_target(MANTA_PYPILE
                  ALL
                  COMMAND ${PYTHON_EXECUTABLE} -m compileall ${CMAKE_CURRENT_BINARY_DIR}/pystage)

include("${MANTA_GLOBALS_CMAKE}") # get MANTA_LIBRARY_PERMISSIONS
set(PYGLOBS "*.py" "*.pyc")
foreach(PYGLOB ${PYGLOBS})
    install (DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pystage" DESTINATION "${INSTALL_TO_DIR}" FILE_PERMISSIONS ${MANTA_LIBRARY_PERMISSIONS} FILES_MATCHING PATTERN ${PYGLOB})
endforeach()

