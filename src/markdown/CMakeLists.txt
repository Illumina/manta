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

function(process_user_guide name)
    set (UG "${name}.md")
    set (UL "${CMAKE_CURRENT_BINARY_DIR}/${name}.html")

    configure_files (${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${UG})

    add_custom_command(
        OUTPUT ${UL}
        COMMAND ${HOEDOWN_PROG} ${CMAKE_CURRENT_BINARY_DIR}/${UG} > ${UL}
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${UG}
        DEPENDS ${THIS_MARKDOWN}
        COMMENT "Converting ${name} to html")

    add_custom_target(${THIS_PROJECT_NAME}_${name} ALL DEPENDS ${UL})

    install(FILES ${UL} DESTINATION ${THIS_DOCDIR}/html)
endfunction()

process_user_guide("mantaUserGuide")

install(FILES "userGuide.css" DESTINATION ${THIS_DOCDIR}/html)
install(FILES "tableOfContents.js" DESTINATION ${THIS_DOCDIR}/html)

