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

#include <array>



template <typename Vec>
void
appendVec(
    Vec& A,
    const Vec& B)
{
    A.insert( A.end(), B.begin(), B.end() );
}



/// an SV candidate with additional details pertaining to input read evidence which is useful for filtration
///
struct FatSVCandidate : public SVCandidate
{
    typedef SVCandidate base_t;

    FatSVCandidate()
        : base_t()
    {}

    FatSVCandidate(const SVCandidate& copy)
        : base_t(copy)
    {}

    FatSVCandidate(const FatSVCandidate&) = default;
    FatSVCandidate& operator=(const FatSVCandidate&) = default;


    bool
    merge(const FatSVCandidate& rhs)
    {
        if (! base_t::merge(rhs)) return false;
        for (unsigned i(0); i<SVEvidenceType::SIZE; ++i)
        {
            appendVec(bp1EvidenceIndex[i],rhs.bp1EvidenceIndex[i]);
            appendVec(bp2EvidenceIndex[i],rhs.bp2EvidenceIndex[i]);
        }
        return true;
    }

#if 0
    bool
    merge(const SVCandidate& rhs)
    {
        if (! base_t::merge(rhs)) return false;

        return true;
    }
#endif

    bool
    evidenceMerge(const FatSVCandidate& rhs)
    {
        if (! base_t::evidenceMerge(rhs)) return false;
        return true;
    }

    void
    clear()
    {
        base_t::clear();
        for (auto& evi : bp1EvidenceIndex) evi.clear();
        for (auto& evi : bp2EvidenceIndex) evi.clear();
    }

    typedef std::array<std::vector<double>,SVEvidenceType::SIZE> evidenceIndex_t;

    evidenceIndex_t bp1EvidenceIndex;
    evidenceIndex_t bp2EvidenceIndex;
};
