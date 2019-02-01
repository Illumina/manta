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

#include "SVCandidateProcessor.hh"

#include <iostream>

#include "blt_util/log.hh"
#include "manta/SVMultiJunctionCandidateUtil.hh"
#include "svgraph/EdgeInfoUtil.hh"
#include "svgraph/SVLocusSet.hh"



SVCandidateProcessor::
SVCandidateProcessor(
    const GSCOptions& opt,
    const SVLocusScanner& readScanner,
    const char* progName,
    const char* progVersion,
    const SVLocusSet& cset,
    std::shared_ptr<EdgeRuntimeTracker> edgeTrackerPtr,
    GSCEdgeStatsManager& edgeStatMan) :
    _opt(opt),
    _cset(cset),
    _edgeTrackerPtr(edgeTrackerPtr),
    _edgeStatMan(edgeStatMan),
    _svRefine(opt, cset.getBamHeader(), cset.getAllSampleReadCounts(), _edgeTrackerPtr),
    _svWriter(opt, readScanner, cset.getBamHeader(), progName, progVersion)
{}



void
SVCandidateProcessor::
evaluateCandidates(
    const EdgeInfo& edge,
    const std::vector<SVMultiJunctionCandidate>& mjSVs,
    const SVCandidateSetData& svData,
    SupportSamples& svSupports)
{
    const bool isIsolatedEdge(testIsolatedEdge(_cset,edge));

    bool isFindLargeInsertions(isIsolatedEdge);
    if (isFindLargeInsertions)
    {
        for (const SVMultiJunctionCandidate& mjCandidateSV : mjSVs)
        {
            for (const SVCandidate& candidateSV : mjCandidateSV.junction)
            {
                if (! isComplexSV(candidateSV)) isFindLargeInsertions=false;
            }
        }
    }

    _svRefine.clearEdgeData();
    for (const auto& cand : mjSVs)
    {
        evaluateCandidate(edge,cand,svData,isFindLargeInsertions, svSupports);
    }
}



void
SVCandidateProcessor::
evaluateCandidate(
    const EdgeInfo& edge,
    const SVMultiJunctionCandidate& mjCandidateSV,
    const SVCandidateSetData& svData,
    const bool isFindLargeInsertions,
    SupportSamples& svSupports)
{
    assert(! mjCandidateSV.junction.empty());

    const unsigned junctionCount(mjCandidateSV.junction.size());

    if (_opt.isVerbose)
    {
        log_os << __FUNCTION__ << ": Starting analysis for SV candidate containing " << junctionCount << " junctions. Low-resolution junction candidate ids:";
        for (const SVCandidate& sv : mjCandidateSV.junction)
        {
            log_os << " " << sv.candidateIndex;
        }
        log_os << "\n";
    }
#ifdef DEBUG_GSV
    log_os << __FUNCTION__ << ": CandidateSV: " << mjCandidateSV << "\n";
#endif


    const bool isComplex(isComplexSV(mjCandidateSV));
    _edgeTrackerPtr->addCandidate(isComplex);

    _edgeStatMan.updateJunctionCandidateCounts(edge, junctionCount, isComplex);

    // true if any junction returns a complex (indel) assembly result
    //
    // if true, and if the candidate is multi-junction then score and write each junction independently:
    bool isAnySmallAssembler(false);
    std::vector<SVCandidateAssemblyData> mjAssemblyData(junctionCount);

    if (! _opt.isSkipAssembly)
    {
        const TimeScoper assmTime(_edgeTrackerPtr->assemblyTime);
        for (unsigned junctionIndex(0); junctionIndex<junctionCount; ++junctionIndex)
        {
            const SVCandidate& candidateSV(mjCandidateSV.junction[junctionIndex]);
            SVCandidateAssemblyData& assemblyData(mjAssemblyData[junctionIndex]);
            try
            {
                _svRefine.getCandidateAssemblyData(candidateSV, isFindLargeInsertions, assemblyData);
            }
            catch (illumina::common::ExceptionData& e)
            {
                std::ostringstream oss;
                oss << "Exception caught while attempting to assemble " << candidateSV;
                e << boost::error_info<struct assembly_candidate_info,std::string>(oss.str());
                throw;
            }
            catch (...)
            {
                log_os << "Exception caught while attempting to assemble " << candidateSV << "\n";
                throw;
            }


            if (_opt.isVerbose)
            {
                log_os << __FUNCTION__ << ": Candidate assembly complete for junction " << junctionIndex << "/" << junctionCount << ". Assembled candidate count: " << assemblyData.svs.size() << "\n";
            }

            if (! assemblyData.svs.empty())
            {
                const unsigned assemblyCount(assemblyData.svs.size());

                const bool isSpanning(assemblyData.isSpanning);

                _edgeStatMan.updateAssemblyCount(edge, assemblyCount, isSpanning);

                if (isSpanning)
                {
                    // can't be multi-junction and multi-assembly at the same time:
                    assert(! ((junctionCount>1) && (assemblyCount>1)));
                }
                else
                {
                    isAnySmallAssembler=true;
                }

                // fill in assembly tracking data:
                _edgeTrackerPtr->addAssembledCandidate(isComplex);
            }
            else
            {
                _edgeStatMan.updateAssemblyCount(edge, 0, assemblyData.isSpanning, assemblyData.isOverlapSkip);
            }
        }
    }

    SVMultiJunctionCandidate mjAssembledCandidateSV;
    mjAssembledCandidateSV.junction.resize(junctionCount);
    std::vector<bool> isJunctionFiltered(junctionCount,false);

    std::vector<unsigned> junctionTracker(junctionCount,0);
    while (true)
    {
        // Note this loop is an accident -- it was intended to enumerate all assembly
        // combinations for  multiple junctions with multiple assemblies each.
        // It doesn't do that -- but the broken thing it does, in fact, do, is what we
        // want for the isAnySmallAssembler case so it's well enough for now.
        //
        /// TODO: update docs -- what is it we "want for the isAnySmallAssembler case"?
        bool isWrite(false);
        for (unsigned junctionIndex(0); junctionIndex<junctionCount; ++junctionIndex)
        {
            const SVCandidateAssemblyData& assemblyData(mjAssemblyData[junctionIndex]);
            unsigned& assemblyIndex(junctionTracker[junctionIndex]);

            if (assemblyData.svs.empty())
            {
                if (assemblyIndex != 0) continue;
#ifdef DEBUG_GSV
                log_os << __FUNCTION__ << ": score and output low-res candidate junction " << junctionIndex << "\n";
#endif
                mjAssembledCandidateSV.junction[junctionIndex] = mjCandidateSV.junction[junctionIndex];
            }
            else
            {
                if (assemblyIndex >= assemblyData.svs.size()) continue;
                const SVCandidate& assembledSV(assemblyData.svs[assemblyIndex]);
#ifdef DEBUG_GSV
                log_os << __FUNCTION__ << ": score and output assembly candidate junction " << junctionIndex << ": " << assembledSV << "\n";
#endif
                mjAssembledCandidateSV.junction[junctionIndex] = assembledSV;
            }
            assemblyIndex++;
            isWrite = true;
        }
        if (! isWrite) break;


        // if any junctions go into the small assembler (for instance b/c the breakends are too close), then
        // treat all junctions as independent:
        //
        const TimeScoper scoreTime(_edgeTrackerPtr->scoreTime);
        if ((junctionCount>1) && isAnySmallAssembler)
        {
            // call each junction independently by providing a filter vector in each iteration
            // which only leaves a single junction unfiltered:
            for (unsigned junctionIndex(0); junctionIndex<junctionCount; ++junctionIndex)
            {
                std::vector<bool> isJunctionFilteredHack(junctionCount,true);
                isJunctionFilteredHack[junctionIndex] = false;
                _svWriter.writeSV(edge, svData, mjAssemblyData, mjAssembledCandidateSV,
                                  isJunctionFilteredHack, svSupports);
            }
        }
        else
        {
            _svWriter.writeSV(edge, svData, mjAssemblyData, mjAssembledCandidateSV,
                              isJunctionFiltered, svSupports);
        }
    }
}
