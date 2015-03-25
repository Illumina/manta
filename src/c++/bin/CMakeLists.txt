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
## Configuration file for the c++/bin subdirectory
##
## author Come Raczy
##
################################################################################

include(${THIS_CXX_EXECUTABLE_CMAKE})

file (GLOB THIS_PROGRAM_SOURCE_LIST [a-zA-Z0-9]*.cpp)

##
## Generic rule for all the other programs
##
foreach(THIS_PROGRAM_SOURCE ${THIS_PROGRAM_SOURCE_LIST})
    get_filename_component(THIS_PROGRAM ${THIS_PROGRAM_SOURCE} NAME_WE)
    set(THIS_APPLICATION_LIB manta_${THIS_PROGRAM})
    add_executable        (${THIS_PROGRAM} ${THIS_PROGRAM_SOURCE})
    target_link_libraries (${THIS_PROGRAM}  ${THIS_APPLICATION_LIB} ${THIS_AVAILABLE_LIBRARIES}
                           ${SAMTOOLS_LIBRARY} ${HTSLIB_LIBRARY} ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES}
                           ${THIS_ADDITIONAL_LIB})
    install(TARGETS ${THIS_PROGRAM} RUNTIME DESTINATION ${THIS_LIBEXECDIR})
endforeach()

