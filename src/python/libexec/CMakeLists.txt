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

set(INSTALL_TO_DIR "${THIS_LIBEXECDIR}")

file(RELATIVE_PATH THIS_RELATIVE_PYTHON_LIBDIR "${INSTALL_TO_DIR}" "${THIS_PYTHON_LIBDIR}")
file(RELATIVE_PATH THIS_RELATIVE_LIBEXECDIR "${INSTALL_TO_DIR}" "${THIS_LIBEXECDIR}")

include("${THIS_MACROS_CMAKE}")
configure_files("${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}" "*.py")

include("${THIS_GLOBALS_CMAKE}") # get THIS_EXECUTABLE_PERMISSIONS
install_fileglob("${CMAKE_CURRENT_BINARY_DIR}" "${INSTALL_TO_DIR}" "*.py"
                 "${THIS_EXECUTABLE_PERMISSIONS}")

