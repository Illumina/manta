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

#include "EdgeRuntimeTracker.hh"
#include "GSCEdgeStatsManager.hh"
#include "GSCOptions.hh"
#include "SVCandidateAssemblyRefiner.hh"
#include "SVScorer.hh"

#include "common/OutStream.hh"
#include "manta/JunctionIdGenerator.hh"
#include "manta/SVCandidateAssemblyData.hh"
#include "manta/SVMultiJunctionCandidate.hh"
#include "format/VcfWriterCandidateSV.hh"
#include "format/VcfWriterDiploidSV.hh"
#include "format/VcfWriterSomaticSV.hh"
#include "truth/TruthTracker.hh"

#include <memory>

//#define DEBUG_GSV



struct SVWriter
{
    SVWriter(
        const GSCOptions& initOpt,
        const SVLocusSet& cset,
        const char* progName,
        const char* progVersion,
        TruthTracker& truthTracker);

    void
    writeSV(
        const EdgeInfo& edge,
        const SVCandidateSetData& svData,
        const std::vector<SVCandidateAssemblyData>& assemblyData,
        const SVMultiJunctionCandidate& mjSV,
        const std::vector<bool>& isInputJunctionFiltered);

    ///////////////////////// data:
    const GSCOptions& opt;
    const bool isSomatic;

    SVScorer svScore;
    std::vector<SVModelScoreInfo> mjModelScoreInfo;

    OutStream candfs;
    OutStream dipfs;
    OutStream somfs;

    VcfWriterCandidateSV candWriter;
    VcfWriterDiploidSV diploidWriter;
    VcfWriterSomaticSV somWriter;

    TruthTracker& _truthTracker;

    JunctionIdGenerator _idgen;
};


struct SVCandidateProcessor
{
    SVCandidateProcessor(
        const GSCOptions& opt,
        const char* progName,
        const char* progVersion,
        const SVLocusSet& cset,
        TruthTracker& truthTracker,
        EdgeRuntimeTracker& edgeTracker,
        GSCEdgeStatsManager& _edgeStatMan);

    void
    evaluateCandidates(
        const EdgeInfo& edge,
        const std::vector<SVMultiJunctionCandidate>& mjSVs,
        const SVCandidateSetData& svData);

private:

    void
    evaluateCandidate(
        const EdgeInfo& edge,
        const SVMultiJunctionCandidate& mjCandidateSV,
        const SVCandidateSetData& svData,
        const bool isFindLargeInsertions);

    const GSCOptions& _opt;
    const SVLocusSet& _cset;
    TruthTracker& _truthTracker;
    EdgeRuntimeTracker& _edgeTracker;
    GSCEdgeStatsManager& _edgeStatMan;
    SVCandidateAssemblyRefiner _svRefine;
    SVWriter _svWriter;
};
