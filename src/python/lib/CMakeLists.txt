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

include ("${MANTA_MACROS_CMAKE}")
configure_files ("${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}" "*.py")

add_custom_target(MANTA_PYPILE
                  ALL
                  COMMAND ${PYTHON_EXECUTABLE} -m compileall ${CMAKE_CURRENT_BINARY_DIR})

include ("${MANTA_GLOBALS_CMAKE}") # get MANTA_LIBRARY_PERMISSIONS
set (PYGLOBS "*.py" "*.pyc")
foreach (PYGLOB ${PYGLOBS})
    install_files ("${CMAKE_CURRENT_BINARY_DIR}" "${MANTA_PYTHON_LIBDIR}" ${PYGLOB}
                   "${MANTA_LIBRARY_PERMISSIONS}")
endforeach ()

