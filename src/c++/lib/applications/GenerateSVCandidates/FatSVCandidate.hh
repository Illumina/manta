//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Chris Saunders
///

#pragma once

#include "manta/SVCandidate.hh"

#include <array>
#include <iosfwd>

/// extend vector A with the contents of B
///
/// example:
/// vector<int> A = {1,2};
/// vector<int> B = {2,3};
/// appendVec(A,B);
/// assert(A == {1,2,2,3});
///
template <typename Vec>
void
appendVec(
    Vec& A,
    const Vec& B)
{
    A.insert( A.end(), B.begin(), B.end() );
}



/// \brief An SV candidate with additional details pertaining to input read evidence
///
/// The extra read evidence provided in this version of SV candidate is useful for filtration
///
struct FatSVCandidate : public SVCandidate
{
    typedef SVCandidate base_t;

    FatSVCandidate()
        : base_t()
    {}

    explicit
    FatSVCandidate(const SVCandidate& copy)
        : base_t(copy)
    {}

    FatSVCandidate(const FatSVCandidate&) = default;
    FatSVCandidate& operator=(const FatSVCandidate&) = default;


    bool
    merge(
        const FatSVCandidate& rhs,
        const bool isExpandRegion = true)
    {
        if (! base_t::merge(rhs, isExpandRegion)) return false;
        for (unsigned evidenceTypeIndex(0); evidenceTypeIndex<SVEvidenceType::SIZE; ++evidenceTypeIndex)
        {
            appendVec(bp1EvidenceIndex[evidenceTypeIndex],rhs.bp1EvidenceIndex[evidenceTypeIndex]);
            appendVec(bp2EvidenceIndex[evidenceTypeIndex],rhs.bp2EvidenceIndex[evidenceTypeIndex]);
        }
        return true;
    }

    /// a 2d array type to track breakpoint evidence, the first dimension is evidence type
    /// and the inner dimension is a vector with size equal to the number of (confident-mapping) observations.
    /// For each observation the inner-dimension value provides the index of the read used as an observation, which
    /// can be used to estimate signal density vs. all reads.
    typedef std::array<std::vector<double>,SVEvidenceType::SIZE> evidenceIndex_t;

    evidenceIndex_t bp1EvidenceIndex;
    evidenceIndex_t bp2EvidenceIndex;
};


std::ostream&
operator<<(std::ostream& os, const FatSVCandidate& svc);
