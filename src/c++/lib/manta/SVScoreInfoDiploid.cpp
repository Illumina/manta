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

#include "manta/SVScoreInfoDiploid.hh"
#include "blt_util/log.hh"

#include "boost/foreach.hpp"

#include <iostream>



std::ostream&
operator<<(
    std::ostream& os,
    const SVScoreInfoDiploid& sid)
{
    os << "DiploidSVScoreInfo bp1MaxDepth=" << sid.bp1MaxDepth << " bp2MaxDepth=" << sid.bp2MaxDepth
       << " altScore=" << sid.altScore
       << " gtScore=" << sid.gtScore
       << " gt=" << DIPLOID_GT::label(sid.gt)
       << "\n";
    os << "Normal sample info " << sid.normal;
    BOOST_FOREACH(const std::string& filter, sid.filters)
    {
        os << " " << filter;
    }
    return os;
}
