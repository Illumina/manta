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
## CMake configuration file for common installation macros
##
## authors: Roman Petrovski, Mauricio Varea
##
################################################################################

macro(configure_files srcDir destDir pattern)
    file(GLOB templateFiles RELATIVE ${srcDir} ${srcDir}/${pattern})
    foreach(templateFile ${templateFiles})
        #message(STATUS "Configuring file ${srcDir}/${templateFile}")
        configure_file(${srcDir}/${templateFile} ${destDir}/${templateFile} @ONLY)
    endforeach()
endmacro()


macro(install_fileglob srcDir destDir pattern perm)
    file(GLOB templateFiles ${srcDir}/${pattern})
    install(FILES ${templateFiles} DESTINATION ${destDir} PERMISSIONS ${perm})
endmacro()


#
# Macro to find libraries, with support for static-only search
#
macro(static_find_library name header library)
    if    (NOT ${name}_INCLUDE_DIR)
        find_path(${name}_INCLUDE_DIR ${header}
                  HINTS ENV C_INCLUDE_PATH ENV CPATH ENV CPLUS_INCLUDE_PATH)
    endif ()
    if    (${name}_INCLUDE_DIR AND NOT ${name}_LIBRARY)
        find_library(${name}_LIBRARY
                     NAMES "${LIBRARY_PREFIX}${library}${LIBRARY_SUFFIX}"
                     HINTS ENV LIBRARY_PATH)
    endif ()
    if(${name}_INCLUDE_DIR AND ${name}_LIBRARY)
        set (HAVE_${name} ${${name}_LIBRARY})
        message (STATUS "Found ${name}  header: ${${name}_INCLUDE_DIR}/${header}")
        message (STATUS "Found ${name} library: ${${name}_LIBRARY}")
    endif()
endmacro()



#
# get all sub directories: (refaim@stackoverflow)
#
function(subdirlist result curdir)
    file(GLOB children RELATIVE ${curdir} ${curdir}/*)
    set(dirlist "")
    foreach(child ${children})
        if(IS_DIRECTORY ${curdir}/${child})
            set(dirlist ${dirlist} ${child})
        endif()
    endforeach()
    set(${result} ${dirlist} PARENT_SCOPE)
endfunction()
