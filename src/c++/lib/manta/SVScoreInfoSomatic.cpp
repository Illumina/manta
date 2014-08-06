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

#include "manta/SVScoreInfoSomatic.hh"
#include "blt_util/log.hh"

#include <iostream>



std::ostream&
operator<<(
    std::ostream& os,
    const SVScoreInfoSomatic& sis)
{
    os << "SomaticSVScoreInfo somaticScore: " << sis.somaticScore << " sstier: " << sis.somaticScoreTier << " filters: ";
    for (const std::string& filter : sis.filters)
    {
        os << " " << filter;
    }
    return os;
}
