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

set(INSTALL_TO_DIR "${THIS_BINDIR}")

file(RELATIVE_PATH THIS_RELATIVE_PYTHON_LIBDIR "${INSTALL_TO_DIR}" "${THIS_PYTHON_LIBDIR}")

include ("${THIS_MACROS_CMAKE}")

configure_files("${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}" "*.py")
configure_files("${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}" "*.ini")

include("${THIS_GLOBALS_CMAKE}") # required for *_PERMISSIONS values

install_fileglob("${CMAKE_CURRENT_BINARY_DIR}" "${INSTALL_TO_DIR}" "*.py"
                 "${THIS_EXECUTABLE_PERMISSIONS}")
install_fileglob("${CMAKE_CURRENT_BINARY_DIR}" "${INSTALL_TO_DIR}" "*.ini"
                 "${THIS_LIBRARY_PERMISSIONS}")

