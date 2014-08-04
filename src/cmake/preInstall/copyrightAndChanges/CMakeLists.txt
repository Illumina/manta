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

#
# Configuration file for the COPYRIGHT and Changes file installation
#
# author Roman Petrovski
#

# Installing top level components
foreach(topfile ChangeLog.txt COPYRIGHT.txt Illumina_Open_Source_License_1.pdf)
    install(FILES "${CMAKE_SOURCE_DIR}/../${topfile}" DESTINATION "${THIS_DATADIR}")
endforeach()