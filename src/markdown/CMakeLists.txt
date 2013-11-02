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

################################################################################
##
## file CMakeLists.txt
##
## Configuration file for the markdown subfolder
##
## author Roman Petrovski
##
################################################################################
include ("${MANTA_GLOBALS_CMAKE}")

##
## build the documentation when available
##
include (FindMarkdown)
include (UseMarkdown)
message (STATUS "Markdown: ${MARKDOWN_EXECUTABLE}.")
if (MARKDOWN_FOUND)
    set (MUG "mantaUserGuide.md")
    configure_files (${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${MUG})

    set(MARKDOWN_PARAMS -x tables -x toc)
    add_markdown_target(mantaUserGuide ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${MUG})
    install_markdown_target(mantaUserGuide DESTINATION ${MANTA_DOCDIR}/html)
endif ()



