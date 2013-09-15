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
// <https://github.com/sequencing/licenses/>
//

///
/// \author Ole Schulz-Trieglaff
///

#include "AssembledContig.hh"

#include <iostream>



std::ostream&
operator<<(std::ostream& os, const AssembledContig& contig)
{
    os << "CONTIG size: " << contig.seq.size()
       << " seedCount: " << contig.seedReadCount
       << " seq: " << contig.seq
       << "\n";

    return os;
}
