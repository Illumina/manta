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

set(INSTALL_TO_DIR "${MANTA_BINDIR}")

file(RELATIVE_PATH MANTA_RELATIVE_PYTHON_LIBDIR "${INSTALL_TO_DIR}" "${MANTA_PYTHON_LIBDIR}")

include ("${MANTA_MACROS_CMAKE}")

configure_files("${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}" "*.py")
configure_files("${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}" "*.ini")

include("${MANTA_GLOBALS_CMAKE}") # required for *_PERMISSIONS values

install_fileglob("${CMAKE_CURRENT_BINARY_DIR}" "${INSTALL_TO_DIR}" "*.py"
                 "${MANTA_EXECUTABLE_PERMISSIONS}")
install_fileglob("${CMAKE_CURRENT_BINARY_DIR}" "${INSTALL_TO_DIR}" "*.ini"
                 "${MANTA_LIBRARY_PERMISSIONS}")

