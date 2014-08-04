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

subdirlist(subdirs ${CMAKE_CURRENT_SOURCE_DIR})

foreach(THIS_LIB_DIR ${subdirs})
    add_subdirectory(${THIS_LIB_DIR})
endforeach()
