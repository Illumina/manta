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

#include "boost/foreach.hpp"

#include <iostream>



std::ostream&
operator<<(
    std::ostream& os,
    const SVScoreInfoSomatic& sis)
{
    os << "SomaticSVScoreInfo bp1MaxDepth=" << sis.bp1MaxDepth << " bp2MaxDepth=" << sis.bp2MaxDepth << " somaticScore=" << sis.somaticScore << "\n";
    os << "Tumor sample info " << sis.tumor;
    os << "Normal sample info " << sis.normal;
    BOOST_FOREACH(const std::string& filter, sis.filters)
    {
        os << " " << filter;
    }
    return os;
}
