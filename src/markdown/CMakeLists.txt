#
# Manta
# Copyright (c) 2013-2015 Illumina, Inc.
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

    # build-time configuration:
    set (UG_TARGET "${THIS_PROJECT_NAME}_${name}_config")

    add_custom_target(
        ${UG_TARGET}
        DEPENDS ${THIS_BUILDTIME_CONFIG_TARGET}
        COMMAND ${CMAKE_COMMAND}
        -D CONFIG_FILE=${THIS_BUILDTIME_CONFIG_FILE}
        -D SOURCE_FILE=${CMAKE_CURRENT_SOURCE_DIR}/${UG}
        -D DEST_FILE=${CMAKE_CURRENT_BINARY_DIR}/${UG}
        -P ${THIS_MODULE_DIR}/buildTimeConfigure.cmake)

    add_custom_command(
        OUTPUT ${UL}
        COMMAND ${HOEDOWN_PROG} ${CMAKE_CURRENT_BINARY_DIR}/${UG} > ${UL}
        DEPENDS ${UG_TARGET}
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${UG}
        DEPENDS ${THIS_MARKDOWN}
        COMMENT "Converting ${name} to html")

    add_custom_target(${THIS_PROJECT_NAME}_${name} ALL DEPENDS ${UL})

    install(FILES ${UL} DESTINATION ${THIS_DOCDIR}/html)
endfunction()

process_user_guide("mantaUserGuide")

install(FILES "userGuide.css" DESTINATION ${THIS_DOCDIR}/html)
install(FILES "tableOfContents.js" DESTINATION ${THIS_DOCDIR}/html)

