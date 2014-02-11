#!/usr/bin/env bash
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
# get version number from git describe
#

set -o nounset

script_dir=$(dirname $0)

version="UNKNOWN"
cd $script_dir
git_version=$(git describe | sed "s/^v//" 2> /dev/null)
if [ $? == 0 ]; then
    version=$git_version
fi

echo -n $version

