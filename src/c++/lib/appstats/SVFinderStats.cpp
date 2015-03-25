// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
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

#include "SVFinderStats.hh"
#include <iostream>



void
SVFinderStats::
report(std::ostream& os) const
{
    os << "EdgeFilter\t" << edgeFilter << "\n";
    os << "SemiMappedFilter\t" << semiMappedFilter << "\n";
    os << "ComplexLowCountFilter\t" << ComplexLowCountFilter << "\n";
    os << "ComplexLowSignalFilter\t" << ComplexLowSignalFilter << "\n";
}
