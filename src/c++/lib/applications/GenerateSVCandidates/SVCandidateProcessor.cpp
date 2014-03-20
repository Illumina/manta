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

#include "SVCandidateProcessor.hh"
#include "manta/SVMultiJunctionCandidateUtil.hh"

#include "blt_util/log.hh"

#include "boost/foreach.hpp"

#include <iostream>


//#define DEBUG_GSV



SVWriter::
SVWriter(
    const GSCOptions& initOpt,
    const SVLocusSet& cset,
    const char* progName,
    const char* progVersion,
    TruthTracker& truthTracker) :
    opt(initOpt),
    isSomatic(! opt.somaticOutputFilename.empty()),
    svScore(opt, cset.header),
    candfs(opt.candidateOutputFilename),
    dipfs(opt.diploidOutputFilename),
    somfs(opt.somaticOutputFilename),
    candWriter(opt.referenceFilename,cset,candfs.getStream()),
    diploidWriter(opt.diploidOpt, (! opt.chromDepthFilename.empty()),
                  opt.referenceFilename,cset,dipfs.getStream()),
    somWriter(opt.somaticOpt, (! opt.chromDepthFilename.empty()),
              opt.referenceFilename,cset,somfs.getStream()),
    _truthTracker(truthTracker)
{
    if (0 == opt.edgeOpt.binIndex)
    {
        candWriter.writeHeader(progName, progVersion);
        diploidWriter.writeHeader(progName, progVersion);
        if (isSomatic) somWriter.writeHeader(progName, progVersion);
    }
}



static
bool
isAnyFalse(
    const std::vector<bool>& vb)
{
    BOOST_FOREACH(const bool val, vb)
    {
        if (! val) return true;
    }
    return false;
}



void
SVWriter::
writeSV(
    const EdgeInfo& edge,
    const SVCandidateSetData& svData,
    const std::vector<SVCandidateAssemblyData>& mjAssemblyData,
    const SVMultiJunctionCandidate& mjSV,
    const std::vector<bool>& isInputJunctionFiltered)
{
    const unsigned junctionCount(mjSV.junction.size());

    // track filtration for each junction:
    std::vector<bool> isJunctionFiltered(isInputJunctionFiltered);

    // early SV filtering:
    //
    // 2 junction filter types:
    // 1) tests where the junction can fail independently
    // 2) tests where all junctions have to fail for the candidate to be filtered:

    bool isCandidateSpanFail(true);

    for (unsigned junctionIndex(0); junctionIndex<junctionCount; ++junctionIndex)
    {
        const SVCandidateAssemblyData& assemblyData(mjAssemblyData[junctionIndex]);
        const SVCandidate& sv(mjSV.junction[junctionIndex]);

        const bool isCandidateSpanning(assemblyData.isCandidateSpanning);

#ifdef DEBUG_GSV
        log_os << __FUNCTION__ << ": isSpanningSV junction: " <<  isCandidateSpanning << "\n";
#endif

        // junction dependent tests:
        //
        bool isJunctionSpanFail(false);
        if (isCandidateSpanning)
        {
            static const unsigned minCandidateSpanningCount(3);
            if (sv.bp1.getSpanningCount() < minCandidateSpanningCount)
            {
                isJunctionSpanFail=true;
            }
        }
        if (! isJunctionSpanFail) isCandidateSpanFail=false;

        // independent tests -- as soon as one of these fails, we can continue:
        //
        if (! isCandidateSpanning)
        {
            if (sv.isImprecise())
            {
                // in this case a non-spanning low-res candidate went into assembly but
                // did not produce a successful contig alignment:
#ifdef DEBUG_GSV
                log_os << __FUNCTION__ << ": Rejecting candidate junction: imprecise non-spanning SV\n";
#endif
                _truthTracker.reportOutcome(SVLog::IMPRECISE_NON_SPANNING);
                isJunctionFiltered[junctionIndex] = true;
                continue;
            }
        }

        // check min size for candidate output:
        if (isSVBelowMinSize(sv,opt.scanOpt.minCandidateVariantSize))
        {
#ifdef DEBUG_GSV
            log_os << __FUNCTION__ << ": Filtering out candidate below min size before candidate output stage\n";
#endif
            return;
        }
    }

    // revisit dependent tests:
    //
    if (isCandidateSpanFail)
    {
        for (unsigned junctionIndex(0); junctionIndex<junctionCount; ++junctionIndex)
        {
#ifdef DEBUG_GSV
            log_os << __FUNCTION__ << ": Rejecting candidate junction: minCandidateSpanningCount\n";
#endif
            _truthTracker.reportOutcome(SVLog::LOW_SPANNING_COUNT);
            isJunctionFiltered[junctionIndex] = true;
        }
    }

    // check to see if all junctions are filtered, if so skip the whole candidate:
    //
    if (! isAnyFalse(isJunctionFiltered))
    {
#ifdef DEBUG_GSV
        log_os << __FUNCTION__ << ": Rejecting candidate, all junctions filtered.\n";
#endif
        return;
    }

    std::vector<SVId> junctionSVId(junctionCount);

    // write out candidates for each junction independently:
    for (unsigned junctionIndex(0); junctionIndex<junctionCount; ++junctionIndex)
    {
        if (isJunctionFiltered[junctionIndex]) continue;

        const SVCandidateAssemblyData& assemblyData(mjAssemblyData[junctionIndex]);
        const SVCandidate& sv(mjSV.junction[junctionIndex]);
        SVId& svId(junctionSVId[junctionIndex]);

        _idgen.getId(edge, sv, svId);

        candWriter.writeSV(svData, assemblyData, sv, svId);
    }

    // check min size for scoring:
    for (unsigned junctionIndex(0); junctionIndex<junctionCount; ++junctionIndex)
    {
        if (isJunctionFiltered[junctionIndex]) continue;

        const SVCandidate& sv(mjSV.junction[junctionIndex]);
        if (isSVBelowMinSize(sv, opt.minScoredVariantSize))
        {
#ifdef DEBUG_GSV
            log_os << __FUNCTION__ << ": Filtering out candidate junction below min size at scoring stage\n";
#endif
            isJunctionFiltered[junctionIndex] = true;
        }
    }

    // check to see if all junctions are filtered before scoring:
    //
    if (! isAnyFalse(isJunctionFiltered)) return;

    bool isMJEvent(false);
    SVModelScoreInfo mjJointModelScoreInfo;
    svScore.scoreSV(svData, mjAssemblyData, mjSV, isJunctionFiltered, isSomatic, mjModelScoreInfo, mjJointModelScoreInfo, isMJEvent);

    unsigned unfilteredJunctionCount(0);
    for (unsigned junctionIndex(0); junctionIndex<junctionCount; ++junctionIndex)
    {
        if (isJunctionFiltered[junctionIndex]) continue;
        unfilteredJunctionCount++;
    }

    // final scored output is treated (mostly) independently for each junction:
    //
    EventInfo event;
    event.junctionCount=unfilteredJunctionCount;

    for (unsigned junctionIndex(0); junctionIndex<junctionCount; ++junctionIndex)
    {
        if (isJunctionFiltered[junctionIndex]) continue;

        const SVCandidateAssemblyData& assemblyData(mjAssemblyData[junctionIndex]);
        const SVCandidate& sv(mjSV.junction[junctionIndex]);
        const SVModelScoreInfo& modelScoreInfo(mjModelScoreInfo[junctionIndex]);

        const SVId& svId(junctionSVId[junctionIndex]);
        if (isMJEvent && event.label.empty())
        {
            event.label = svId.localId;
        }

        const SVScoreInfo& baseInfo(modelScoreInfo.base);
        const SVModelScoreInfo& scoreInfo(isMJEvent ? mjJointModelScoreInfo : modelScoreInfo);
        const SVScoreInfoDiploid& diploidInfo(scoreInfo.diploid);
        const SVScoreInfoSomatic& somaticInfo(scoreInfo.somatic);

        if (diploidInfo.altScore >= opt.diploidOpt.minOutputAltScore || opt.isRNA) /// TODO remove after adding RNA scoring
        {
            diploidWriter.writeSV(svData, assemblyData, sv, svId, baseInfo, diploidInfo, event);
        }

        if (isSomatic)
        {
            if (somaticInfo.somaticScore > opt.somaticOpt.minOutputSomaticScore)
            {
                somWriter.writeSV(svData, assemblyData, sv, svId, baseInfo, somaticInfo, event);
                _truthTracker.reportOutcome(SVLog::WRITTEN);
            }
            else
            {
                _truthTracker.reportOutcome(SVLog::LOW_SOMATIC_SCORE);
            }
        }
    }
}



SVCandidateProcessor::
SVCandidateProcessor(
    const GSCOptions& opt,
    const char* progName,
    const char* progVersion,
    const SVLocusSet& cset,
    TruthTracker& truthTracker,
    EdgeRuntimeTracker& edgeTracker) :
    _opt(opt),
    _truthTracker(truthTracker),
    _edgeTracker(edgeTracker),
    _svRefine(opt, cset.header),
    _svWriter(opt, cset, progName, progVersion, truthTracker)
{}



void
SVCandidateProcessor::
evaluateCandidate(
    const EdgeInfo& edge,
    const SVMultiJunctionCandidate& mjCandidateSV,
    const SVCandidateSetData& svData)
{
    assert(! mjCandidateSV.junction.empty());

    const unsigned junctionCount(mjCandidateSV.junction.size());

    if (_opt.isVerbose)
    {
        log_os << __FUNCTION__ << ": Starting analysis for SV candidate containing " << junctionCount << " junctions. Low-resolution junction candidate ids:";
        BOOST_FOREACH(const SVCandidate& sv, mjCandidateSV.junction)
        {
            log_os << " " << sv.candidateIndex;
        }
        log_os << "\n";
    }
#ifdef DEBUG_GSV
    log_os << __FUNCTION__ << ": CandidateSV: " << mjCandidateSV << "\n";
#endif


    const bool isComplex(isComplexSV(mjCandidateSV));
    _edgeTracker.addCand(isComplex);

    /// assemble each junction independently:
    bool isAnySmallAssembler(false);
    std::vector<SVCandidateAssemblyData> mjAssemblyData(junctionCount);

    if (! _opt.isSkipAssembly)
    {
        for (unsigned junctionIndex(0); junctionIndex<junctionCount; ++junctionIndex)
        {
            const SVCandidate& candidateSV(mjCandidateSV.junction[junctionIndex]);
            SVCandidateAssemblyData& assemblyData(mjAssemblyData[junctionIndex]);
            _svRefine.getCandidateAssemblyData(candidateSV, svData, assemblyData, _opt.isRNA);

            if (_opt.isVerbose)
            {
                log_os << __FUNCTION__ << ": Candidate assembly complete for junction " << junctionIndex << "/" << junctionCount << ". Assembled candidate count: " << assemblyData.svs.size() << "\n";
            }

            if(! assemblyData.svs.empty())
            {
                const unsigned assemblyCount(assemblyData.svs.size());

                const bool isSpanning(assemblyData.isSpanning);

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
                _edgeTracker.addAssm(isComplex);
                _truthTracker.reportNumAssembled(assemblyData.svs.size());
                for (unsigned assemblyIndex(0); assemblyIndex<assemblyCount; ++assemblyIndex)
                {
                    _truthTracker.addAssembledSV();
                }
            }
        }
    }

    SVMultiJunctionCandidate mjAssembledCandidateSV;
    mjAssembledCandidateSV.junction.resize(junctionCount);
    std::vector<bool> isJunctionFiltered(junctionCount,false);

    std::vector<unsigned> junctionTracker(junctionCount,0);
    while(true)
    {
        /// note this loop is stupid -- it was originally written with the intention of
        /// combinatorially enumerating all possible assembly combinations for the case
        /// of multiple junctions with multiple assemblies each.
        /// It doesn't do that -- but the broken thing it does, in fact, do, is what we want for the
        /// isAnySmallAssembler case so it's well enough for now.
        ///
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


        /// if any junctions go into the small assembler (for instance b/c the breakends are too close), then
        /// treat all junctions as independent:
        ///
        if ((junctionCount>1) && isAnySmallAssembler)
        {
            for (unsigned junctionIndex(0); junctionIndex<junctionCount; ++junctionIndex)
            {
                std::vector<bool> isJunctionFilteredHack(junctionCount,true);
                isJunctionFilteredHack[junctionIndex] = false;
                _svWriter.writeSV(edge, svData, mjAssemblyData, mjAssembledCandidateSV, isJunctionFilteredHack);
            }
        }
        else
        {
            _svWriter.writeSV(edge, svData, mjAssemblyData, mjAssembledCandidateSV, isJunctionFiltered);
        }
    }
}
