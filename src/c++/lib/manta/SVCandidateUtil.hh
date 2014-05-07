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

#include "manta/SVCandidate.hh"


/// returns true if sv is below minimum size:
///
bool
isSVBelowMinSize(
    const SVCandidate& sv,
    const unsigned minSize);


namespace SV_TYPE
{
enum index_t
{
    UNKNOWN,
    INTERTRANSLOC,
    INVERSION,
    INDEL,
    TANDUP,
    COMPLEX
};

inline
const char*
label(const index_t idx)
{
    switch (idx)
    {
    case UNKNOWN:
        return "UNKNOWN";
    case INTERTRANSLOC:
        return "INTERTRANSLOC";
    case INVERSION:
        return "INVERSION";
    case INDEL:
        return "INDEL";
    case TANDUP:
        return "TANDUP";
    case COMPLEX:
        return "COMPLEX";
    default:
        return "UNKNOWN";
    }
}

}


SV_TYPE::index_t
getSVType(const SVCandidate& sv);




/// extended SV_TYPE is like SV_TYPE but separates INDEL into INSERT and DELETE states
namespace EXTENDED_SV_TYPE
{

enum index_t
{
    UNKNOWN,
    INTERTRANSLOC,
    INTRATRANSLOC,
    INVERSION,
    INSERT,
    DELETE,
    TANDUP
};

inline
bool
isSVTransloc(const index_t idx)
{
    switch (idx)
    {
    case INTERTRANSLOC:
    case INTRATRANSLOC:
        return true;
    default:
        return false;
    }
}

inline
bool
isSVIndel(const index_t idx)
{
    switch (idx)
    {
    case INSERT:
    case DELETE:
        return true;
    default:
        return false;
    }
}

// provide a shortened label (mostly from the VCF spec)
inline
const char*
label(const index_t idx)
{
    switch (idx)
    {
    case INTERTRANSLOC:
        return "BND";
    case INTRATRANSLOC:
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

EXTENDED_SV_TYPE::index_t
getExtendedSVType(
    const SVCandidate& sv);


inline
bool
isSpanningSV(const SVCandidate& sv)
{
    using namespace SVBreakendState;
    return (isSimpleBreakend(sv.bp1.state) && isSimpleBreakend(sv.bp2.state));
}


/// complex in this case means that we have no specific hypothesis for the SV --
/// it is just a single genomic region for which we schedule local assembly
///
inline
bool
isComplexSV(const SVCandidate& sv)
{
    using namespace SVBreakendState;
    return ((sv.bp1.state == COMPLEX) && (sv.bp2.state == UNKNOWN));
}

