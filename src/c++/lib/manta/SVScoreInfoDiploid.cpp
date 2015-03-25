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

#include "manta/SVScoreInfoDiploid.hh"
#include "blt_util/log.hh"

#include <iostream>



std::ostream&
operator<<(
    std::ostream& os,
    const SVScoreInfoDiploid& sid)
{
    os << "DiploidSVScoreInfo "
       << " altScore=" << sid.altScore
       << " gtScore=" << sid.gtScore
       << " gt=" << DIPLOID_GT::label(sid.gt)
       << "\n";
    for (const std::string& filter : sid.filters)
    {
        os << " " << filter;
    }
    return os;
}
