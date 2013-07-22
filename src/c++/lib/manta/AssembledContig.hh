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

#include <string>
#include <map>

struct AssembledContig {

    AssembledContig() : seedReadCount(0) {}

    // sequence
    std::string seq;
    // reads used for assembly of contig <read_no,mapping position to contig>
    std::map<unsigned,int> contigReads;
    // no of reads containing the seeding kmer
    unsigned seedReadCount;
};

