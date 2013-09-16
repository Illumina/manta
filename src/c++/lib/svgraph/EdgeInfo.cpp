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
/// \author Chris Saunders
///

#include "EdgeInfo.hh"

#include <iostream>



std::ostream&
operator<<(std::ostream& os, const EdgeInfo& ei)
{
    os << "edgeinfo l,n1,n2: " << ei.locusIndex << " " << ei.nodeIndex1 << " " << ei.nodeIndex2 << "\n";
    return os;
}
