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

#include <iosfwd>
#include <string>
#include <vector>
#include <set>


/// \brief data pertaining to a de-novo assembly contig
///
/// stores for each contig the sequence and the number of reads
/// containing its seeding k-mer
///
struct AssembledContig
{
    AssembledContig() : seedReadCount(0) {}

    std::string seq; ///< contigsequence

    // reads used for assembly of contig <read_no,mapping position to contig>
    //std::map<std::string,int> contigReads;

    unsigned seedReadCount; ///< no of reads containing the seeding kmer

    std::set<unsigned> supportReads;
    std::set<unsigned> rejectReads;
};


std::ostream& operator<<(std::ostream& os, const AssembledContig& contig);


typedef std::vector<AssembledContig> Assembly;

