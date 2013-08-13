#
# Manta
# Copyright (c) 2013 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/downloads/sequencing/licenses/>.
#

################################################################################
##
## Configuration file for the c++/bin subdirectory
##
## author Come Raczy
##
################################################################################

include_directories (BEFORE SYSTEM "${SAMTOOLS_DIR}")

include_directories (${MANTA_CXX_ALL_INCLUDES})

include(${MANTA_CXX_EXECUTABLE_CMAKE})

file (GLOB MANTA_PROGRAM_SOURCE_LIST [a-zA-Z0-9]*.cpp)

##
## Generic rule for all the other programs
##
#set(BAM_LIBRARY "${SAMTOOLS_DIR}/libbam.a")
foreach(MANTA_PROGRAM_SOURCE ${MANTA_PROGRAM_SOURCE_LIST})
    get_filename_component(MANTA_PROGRAM ${MANTA_PROGRAM_SOURCE} NAME_WE)
    set(MANTA_APPLICATION_LIB manta_${MANTA_PROGRAM})
    add_executable        (${MANTA_PROGRAM} ${MANTA_PROGRAM_SOURCE})
    target_link_libraries (${MANTA_PROGRAM}  ${MANTA_APPLICATION_LIB} ${MANTA_AVAILABLE_LIBRARIES}
                           ${SAMTOOLS_LIBRARY} ${Boost_LIBRARIES}
                           ${MANTA_ADDITIONAL_LIB})
    install(TARGETS ${MANTA_PROGRAM} RUNTIME DESTINATION ${MANTA_LIBEXECDIR})
endforeach()


