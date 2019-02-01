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
///

#pragma once

#include "EdgeRuntimeTracker.hh"
#include "GSCEdgeStatsManager.hh"
#include "GSCOptions.hh"
#include "SVCandidateAssemblyRefiner.hh"
#include "manta/SVCandidateSetData.hh"
#include "manta/SVLocusScanner.hh"
#include "manta/SVMultiJunctionCandidate.hh"
#include "SVWriter.hh"
#include "SVSupports.hh"



struct SVCandidateProcessor
{
    SVCandidateProcessor(
        const GSCOptions& opt,
        const SVLocusScanner& readScanner,
        const char* progName,
        const char* progVersion,
        const SVLocusSet& cset,
        std::shared_ptr<EdgeRuntimeTracker> edgeTrackerPtr,
        GSCEdgeStatsManager& edgeStatMan);

    /// Refine initial low-resolution candidates using an assembly step, then score and output final SVs
    void
    evaluateCandidates(
        const EdgeInfo& edge,
        const std::vector<SVMultiJunctionCandidate>& mjSVs,
        const SVCandidateSetData& svData,
        SupportSamples& svSupports);

    /// TestSVCandidateProcessor is a friend structure of SVCandidateProcessor. So that it can access private
    /// members of SVCandidateProcessor.
    friend struct TestSVCandidateProcessor;

private:

    void
    evaluateCandidate(
        const EdgeInfo& edge,
        const SVMultiJunctionCandidate& mjCandidateSV,
        const SVCandidateSetData& svData,
        const bool isFindLargeInsertions,
        SupportSamples& svSupports);

    const GSCOptions& _opt;
    const SVLocusSet& _cset;
    std::shared_ptr<EdgeRuntimeTracker> _edgeTrackerPtr;
    GSCEdgeStatsManager& _edgeStatMan;
    SVCandidateAssemblyRefiner _svRefine;
    SVWriter _svWriter;
};
