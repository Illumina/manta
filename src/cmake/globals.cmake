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


install(CODE "

    # With package generator, the location where files are placed is not the location where they will be run.
    # _FULL_ variables are guaranteed valid only at runtime
    set (CPACK_GENERATOR \"${CPACK_GENERATOR}\")
    set (FULL_PREFIX \"\")
    if (NOT CPACK_GENERATOR)
        set (FULL_PREFIX \"$ENV{DESTDIR}\")
    endif()

    get_filename_component(MANTA_FULL_ETCDIR         \"\${FULL_PREFIX}${MANTA_ETCDIR}\" ABSOLUTE)
    get_filename_component(MANTA_FULL_DATADIR        \"\${FULL_PREFIX}${MANTA_DATADIR}\" ABSOLUTE)
    get_filename_component(MANTA_FULL_BINDIR         \"\${FULL_PREFIX}${MANTA_BINDIR}\" ABSOLUTE)
    get_filename_component(MANTA_FULL_LIBDIR         \"\${FULL_PREFIX}${MANTA_LIBDIR}\" ABSOLUTE)
    get_filename_component(MANTA_FULL_LIBEXECDIR     \"\${FULL_PREFIX}${MANTA_LIBEXECDIR}\" ABSOLUTE)
    get_filename_component(MANTA_FULL_PYTHON_LIBDIR  \"\${FULL_PREFIX}${MANTA_PYTHON_LIBDIR}\" ABSOLUTE)

    set(MANTA_VERSION_FULL \"${MANTA_VERSION_FULL}\")
    set(MANTA_VERSION \"${MANTA_VERSION}\")
    set(MANTA_EXECUTABLE_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)
    set(MANTA_LIBRARY_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
    ")

