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
/// \author Chris Saunders
///

#include "EdgeInfo.hh"

#include <iostream>



void
EdgeInfo::
write(std::ostream& os) const
{
    static const char sep(':');
    os << locusIndex << sep << nodeIndex1 << sep << nodeIndex2;
}



std::ostream&
operator<<(std::ostream& os, const EdgeInfo& ei)
{
    os << "edgeinfo locus:node1:node2: ";
    ei.write(os);
    os << '\n';
    return os;
}
