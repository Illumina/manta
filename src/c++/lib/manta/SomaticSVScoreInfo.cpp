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
// <https://github.com/downloads/sequencing/licenses/>.
//

///
/// \author Chris Saunders
///

#include "manta/SomaticSVScoreInfo.hh"

#include "boost/foreach.hpp"

#include <iostream>

std::ostream&
operator<<(std::ostream& os, const SVSampleInfo& si)
{
    os << "SVSampleInfo bp1SpanReads=" << si.bp1SpanReads << " bp2SpanReads=" << si.bp2SpanReads << " spanPairs=" << si.spanPairs << std::endl;
    return os;
}


std::ostream&
operator<<(std::ostream& os, const SomaticSVScoreInfo& ssi)
{
    os << "SomaticSVScoreInfo bp1MaxDepth=" << ssi.bp1MaxDepth << " bp2MaxDepth=" << ssi.bp2MaxDepth << " somaticScore=" << ssi.somaticScore << std::endl;
    os << "Tumor sample info " << ssi.tumor;
    os << "Normal sample info " << ssi.normal;
    BOOST_FOREACH(const std::string& filter, ssi.filters)
    {
        os << " " << filter;
    }
    return os;
}


