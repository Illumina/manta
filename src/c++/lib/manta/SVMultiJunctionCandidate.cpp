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

#include "manta/SVMultiJunctionCandidate.hh"

#include "boost/foreach.hpp"

#include <iostream>



std::ostream&
operator<<(
    std::ostream& os,
    const SVMultiJunctionCandidate& scc)
{
    static const char indent('\t');
    os << "SVComplexCandidate:\n"
       << indent << "total_breakend_junctions: " << scc.junction.size() << "\n";

    BOOST_FOREACH(const SVCandidate& sv, scc.junction)
    {
        os << sv;
    }

    return os;
}
