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

#include "format/VcfWriterCandidateSV.hh"



void
VcfWriterCandidateSV::
writeSV(
    const EdgeInfo& edge,
    const SVCandidateData& svData,
    const std::vector<SVCandidate>& svs)
{
    unsigned svIndex(0);
    BOOST_FOREACH(const SVCandidate& sv, svs)
    {
        if (sv.bp1.pairCount < _minPairCount) continue;
        writeSVCore(edge, svIndex, sv, svData);
        svIndex++;
    }
}
