//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
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
/// \author Naoki Nariai
///

#pragma once

#include <array>
#include <iosfwd>

#include "manta/SVCandidate.hpp"

/// extend vector A with the contents of B
///
/// example:
/// vector<int> A = {1,2};
/// vector<int> B = {2,3};
/// appendVec(A,B);
/// assert(A == {1,2,2,3});
///
template <typename Vec>
void appendVec(Vec& A, const Vec& B)
{
  A.insert(A.end(), B.begin(), B.end());
}

/// \brief An SV candidate with additional details pertaining to input read evidence
///
/// The extra read evidence provided in this version of SV candidate is useful for filtration
///
struct FatSVCandidate : public SVCandidate {
  typedef SVCandidate base_t;

  FatSVCandidate() : base_t() {}

  explicit FatSVCandidate(const SVCandidate& copy, const unsigned bamCount) : base_t(copy)
  {
    for (unsigned evidenceTypeIndex(0); evidenceTypeIndex < SVEvidenceType::SIZE; ++evidenceTypeIndex) {
      bp1EvidenceIndex[evidenceTypeIndex].resize(bamCount);
      bp2EvidenceIndex[evidenceTypeIndex].resize(bamCount);
    }
  }

  FatSVCandidate(const FatSVCandidate&) = default;
  FatSVCandidate& operator=(const FatSVCandidate&) = default;

  bool merge(const FatSVCandidate& rhs, const bool isExpandRegion = true)
  {
    if (!base_t::merge(rhs, isExpandRegion)) return false;
    for (unsigned evidenceTypeIndex(0); evidenceTypeIndex < SVEvidenceType::SIZE; ++evidenceTypeIndex) {
      for (unsigned bamIndex(0); bamIndex < bp1EvidenceIndex[evidenceTypeIndex].size(); ++bamIndex) {
        appendVec(
            bp1EvidenceIndex[evidenceTypeIndex][bamIndex], rhs.bp1EvidenceIndex[evidenceTypeIndex][bamIndex]);
        appendVec(
            bp2EvidenceIndex[evidenceTypeIndex][bamIndex], rhs.bp2EvidenceIndex[evidenceTypeIndex][bamIndex]);
      }
    }
    return true;
  }

  /// a 3d array type to track breakpoint evidence.
  /// The first dimension is evidence type,
  /// the second dimension is bam index with size equal to the number of input bams, and
  /// the third dimension is evidence read with size equal to the number of confident-mapping observations.
  /// For each observation the value provides the index of the read used as an observation,
  /// which can be used to estimate signal density vs. all reads.
  typedef std::array<std::vector<std::vector<double>>, SVEvidenceType::SIZE> evidenceIndex_t;

  evidenceIndex_t bp1EvidenceIndex;
  evidenceIndex_t bp2EvidenceIndex;
};

std::ostream& operator<<(std::ostream& os, const FatSVCandidate& svc);
