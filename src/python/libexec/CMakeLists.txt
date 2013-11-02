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

include ("${MANTA_GLOBALS_CMAKE}") # get MANTA_EXECUTABLE_PERMISSIONS
install_files ("${CMAKE_CURRENT_BINARY_DIR}" "${MANTA_LIBEXECDIR}" "*.py"
               "${MANTA_EXECUTABLE_PERMISSIONS}")

