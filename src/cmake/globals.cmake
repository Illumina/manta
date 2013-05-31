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
## CMake configuration file to identify the configuration of the system
##
## author Roman Petrovski
##
################################################################################

set(MANTA_ORIG_ETCDIR      "${CMAKE_INSTALL_PREFIX}/${MANTA_ETCDIR}")
set(MANTA_ORIG_DATADIR     "${CMAKE_INSTALL_PREFIX}/${MANTA_DATADIR}")
set(MANTA_ORIG_BINDIR      "${CMAKE_INSTALL_PREFIX}/${MANTA_BINDIR}")
set(MANTA_ORIG_LIBDIR      "${CMAKE_INSTALL_PREFIX}/${MANTA_LIBDIR}")
set(MANTA_ORIG_LIBEXECDIR  "${CMAKE_INSTALL_PREFIX}/${MANTA_LIBEXECDIR}")
set(MANTA_ORIG_PERL_LIBDIR "${CMAKE_INSTALL_PREFIX}/${MANTA_PERL_LIBDIR}")
set(MANTA_ORIG_PYTHON_LIBDIR "${CMAKE_INSTALL_PREFIX}/${MANTA_PYTHON_LIBDIR}")

install(CODE "

    # With package generator, the location where files are placed is not the location where they will be run.
    # _FULL_ variables are guaranteed valid only at runtime
    set (CPACK_GENERATOR \"${CPACK_GENERATOR}\")
    if (CPACK_GENERATOR)
        get_filename_component(MANTA_FULL_ETCDIR       \"${MANTA_ORIG_ETCDIR}\" ABSOLUTE)
        get_filename_component(MANTA_FULL_DATADIR      \"${MANTA_ORIG_DATADIR}\" ABSOLUTE)
        get_filename_component(MANTA_FULL_BINDIR       \"${MANTA_ORIG_BINDIR}\" ABSOLUTE)
        get_filename_component(MANTA_FULL_LIBDIR       \"${MANTA_ORIG_LIBDIR}\" ABSOLUTE)
        get_filename_component(MANTA_FULL_LIBEXECDIR   \"${MANTA_ORIG_LIBEXECDIR}\" ABSOLUTE)
        get_filename_component(MANTA_FULL_PERL_LIBDIR  \"${MANTA_ORIG_PERL_LIBDIR}\" ABSOLUTE)
        get_filename_component(MANTA_FULL_PYTHON_LIBDIR  \"${MANTA_ORIG_PYTHON_LIBDIR}\" ABSOLUTE)
    else ()
        get_filename_component(MANTA_FULL_ETCDIR       \"\$ENV{DESTDIR}${MANTA_ORIG_ETCDIR}\" ABSOLUTE)
        get_filename_component(MANTA_FULL_DATADIR      \"\$ENV{DESTDIR}${MANTA_ORIG_DATADIR}\" ABSOLUTE)
        get_filename_component(MANTA_FULL_BINDIR       \"\$ENV{DESTDIR}${MANTA_ORIG_BINDIR}\" ABSOLUTE)
        get_filename_component(MANTA_FULL_LIBDIR       \"\$ENV{DESTDIR}${MANTA_ORIG_LIBDIR}\" ABSOLUTE)
        get_filename_component(MANTA_FULL_LIBEXECDIR   \"\$ENV{DESTDIR}${MANTA_ORIG_LIBEXECDIR}\" ABSOLUTE)
        get_filename_component(MANTA_FULL_PERL_LIBDIR  \"\$ENV{DESTDIR}${MANTA_ORIG_PERL_LIBDIR}\" ABSOLUTE)
        get_filename_component(MANTA_FULL_PYTHON_LIBDIR  \"\$ENV{DESTDIR}${MANTA_ORIG_PYTHON_LIBDIR}\" ABSOLUTE)
    endif ()

    # _DEST_ variables always point to location where files are copied
    get_filename_component(MANTA_DEST_ETCDIR       \"\$ENV{DESTDIR}${MANTA_ORIG_ETCDIR}\" ABSOLUTE)
    get_filename_component(MANTA_DEST_DATADIR      \"\$ENV{DESTDIR}${MANTA_ORIG_DATADIR}\" ABSOLUTE)
    get_filename_component(MANTA_DEST_BINDIR       \"\$ENV{DESTDIR}${MANTA_ORIG_BINDIR}\" ABSOLUTE)
    get_filename_component(MANTA_DEST_LIBDIR       \"\$ENV{DESTDIR}${MANTA_ORIG_LIBDIR}\" ABSOLUTE)
    get_filename_component(MANTA_DEST_LIBEXECDIR   \"\$ENV{DESTDIR}${MANTA_ORIG_LIBEXECDIR}\" ABSOLUTE)
    get_filename_component(MANTA_DEST_PERL_LIBDIR  \"\$ENV{DESTDIR}${MANTA_ORIG_PERL_LIBDIR}\" ABSOLUTE)
    get_filename_component(MANTA_DEST_PYTHON_LIBDIR  \"\$ENV{DESTDIR}${MANTA_ORIG_PYTHON_LIBDIR}\" ABSOLUTE)

    
    set(MANTA_VERSION_MAJOR \"${MANTA_VERSION_MAJOR}\")
    set(MANTA_VERSION_MINOR \"${MANTA_VERSION_MINOR}\")
    set(MANTA_VERSION_PATCH \"${MANTA_VERSION_PATCH}\")
    set(MANTA_VERSION_FULL \"${MANTA_VERSION_FULL}\")
    set(MANTA_EXECUTABLE_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)
    set(MANTA_LIBRARY_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
    ")

