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

#pragma once

#include "svgraph/EdgeInfo.hh"
#include "manta/SVCandidateUtil.hh"
#include "boost/format.hpp"

#include <string>


/// A pair of ids for both ends of a single SV junction
///
/// the mateid will only be defined for tranlocations, and empty otherwise
///
struct SVId
{
    const char*
    getLabel() const
    {
        return EXTENDED_SV_TYPE::label(svType);
    }

    EXTENDED_SV_TYPE::index_t svType = EXTENDED_SV_TYPE::UNKNOWN;
    std::string localId;
    std::string mateId;
};


/// create IDs for each variant that are guaranteed to be unique for a single
/// manta run
///
struct JunctionIdGenerator
{
    JunctionIdGenerator() :
        _SVIdFormatter("Manta%s:%i:%i:%i:%i:%i:%i")
    {}

    void
    getId(
        const EdgeInfo& edge,
        const SVCandidate& sv,
        const bool isRNA,
        SVId& svId);

private:
    boost::format _SVIdFormatter;
};
