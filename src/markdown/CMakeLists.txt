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
    set(MARKDOWN_PARAMS -x tables -x toc)
    add_markdown_target(manual ${CMAKE_CURRENT_BINARY_DIR} "mantaUserGuide.md")
    install_markdown_target(manual DESTINATION ${MANTA_DOCDIR}/html)
endif ()



