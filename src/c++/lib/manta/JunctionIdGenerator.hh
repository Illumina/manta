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

#pragma once

#include "svgraph/EdgeInfo.hh"
#include "manta/SVCandidate.hh"
#include "boost/format.hpp"

#include <string>



/// extended SV_TYPE is like SV_TYPE but separates INDEL into INSERT and DELETE states
namespace EXTENDED_SV_TYPE
{

enum index_t {
    UNKNOWN,
    INTERTRANSLOC,
    INVERSION,
    INSERT,
    DELETE,
    TANDUP
};

// provide a shortened label (mostly from the VCF spec)
inline
const char*
label(const index_t idx)
{
    switch (idx)
    {
    case INTERTRANSLOC:
        return "BND";
    case INVERSION:
        return "INV";
    case INSERT:
        return "INS";
    case DELETE:
        return "DEL";
    case TANDUP:
        return "DUP:TANDEM";
    default:
        return "UNKNOWN";
    }
}
}

/// A pair of ids for both ends of a single SV junction
///
/// the mateid will only be defined for tranlocations, and empty otherwise
///
struct SVId
{
    SVId() :
        svType(EXTENDED_SV_TYPE::UNKNOWN)
    {}

    const char*
    getLabel() const
    {
        return EXTENDED_SV_TYPE::label(svType);
    }

    EXTENDED_SV_TYPE::index_t svType;
    std::string localId;
    std::string mateId;
};


/// create IDs for each variant that are guaranteed to be unique for a single
/// manta run
///
struct JunctionIdGenerator
{
    JunctionIdGenerator() :
        _transLocIdFormatter("MantaBND:%i:%i:%i:%i:"),
        _otherSVIdFormatter("Manta%s:%i:%i:%i:%i:%i:%i")
    {}

    void
    getId(
        const EdgeInfo& edge,
        const SVCandidate& sv,
        SVId& svId);

private:
    boost::format _transLocIdFormatter;
    boost::format _otherSVIdFormatter;
};

