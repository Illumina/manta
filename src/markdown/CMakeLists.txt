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

##
## build the documentation when available
##
set (MUP "mantaUserGuide")
set (MUG "${MUP}.md")
set (MUL "${CMAKE_CURRENT_BINARY_DIR}/${MUP}.html")

configure_files (${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${MUG})

add_custom_TARGET(MANTA_USERGUIDE
                  ALL
                  COMMAND ${HOEDOWN_PROG} ${CMAKE_CURRENT_BINARY_DIR}/${MUG} > ${MUL}
                  DEPENDS MANTA_HOEDOWN
                  COMMENT "Converting user guide to html")

install(FILES ${MUL} DESTINATION ${MANTA_DOCDIR}/html)

