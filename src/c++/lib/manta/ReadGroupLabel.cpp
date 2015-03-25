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

#include "manta/ReadGroupLabel.hh"

#include <iostream>



std::ostream&
operator<<(
    std::ostream& os,
    const ReadGroupLabel& rgl)
{
    os << "read group '" << rgl.rgLabel << "' in bam file '" << rgl.bamLabel;
    return os;
}
