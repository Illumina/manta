// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

///
/// \author Ole Schulz-Trieglaff
///

#pragma once

#include <map>
#include <string>


struct AssemblyRead
{
    AssemblyRead()
        : used(false) {}

    AssemblyRead(std::string s,
                 bool u)
        : seq(s), used(u) {}

    std::string seq;
    bool used;
};


typedef std::map<std::string,AssemblyRead> AssemblyReadMap;
