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
operator<<(std::ostream& os, const SomaticSVScoreInfo& ssi)
{
    BOOST_FOREACH(const std::string& filter, ssi.filters)
    {
        os << " " << filter;
    }
    return os;
}


