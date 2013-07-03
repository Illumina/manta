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

subdirlist(subdirs ${CMAKE_CURRENT_SOURCE_DIR})

foreach(MANTA_LIB_DIR ${subdirs})
    add_subdirectory(${MANTA_LIB_DIR})
endforeach()
