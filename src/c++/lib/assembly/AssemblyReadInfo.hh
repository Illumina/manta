// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Ole Schulz-Trieglaff
///

#pragma once

#include <vector>
#include <string>


/// information added to each read in the process of assembly
///
struct AssemblyReadInfo
{
    AssemblyReadInfo() :
        isUsed(false),
        isFiltered(false),
        contigId(0)
    {}

    bool isUsed;
    bool isFiltered; ///< if true, the read was 'used' but filtered out, so there is no meaningful contig id association
    unsigned contigId; ///< index of the contig that this read is used in
};


typedef std::vector<std::string> AssemblyReadInput;
typedef std::vector<bool> AssemblyReadReversal;
typedef std::vector<AssemblyReadInfo> AssemblyReadOutput;
