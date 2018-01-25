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
/// \author Naoki Nariai
///

#include "SVCandidateProcessor.hh"
#include "manta/SVMultiJunctionCandidateUtil.hh"
#include "svgraph/EdgeInfoUtil.hh"

#include "blt_util/log.hh"

#include <iostream>


//#define DEBUG_GSV



SVWriter::
SVWriter(
    const GSCOptions& initOpt,
    const SVLocusScanner& readScanner,
    const SVLocusSet& cset,
    const char* progName,
    const char* progVersion) :
    opt(initOpt),
    isSomatic(! opt.somaticOutputFilename.empty()),
    isTumorOnly(! opt.tumorOutputFilename.empty()),
    svScore(opt, readScanner, cset.header),
    candfs(opt.candidateOutputFilename),
    dipfs(opt.diploidOutputFilename),
    somfs(opt.somaticOutputFilename),
    tumfs(opt.tumorOutputFilename),
    rnafs(opt.rnaOutputFilename),
    candWriter(opt.referenceFilename, cset, candfs.getStream(),
               opt.isOutputContig),
    diploidWriter(opt.diploidOpt, (! opt.chromDepthFilename.empty()),
                  opt.referenceFilename, cset, dipfs.getStream(),
                  opt.isOutputContig),
    somWriter(opt.somaticOpt, (! opt.chromDepthFilename.empty()),
              opt.referenceFilename, cset,somfs.getStream(),
              opt.isOutputContig),
    tumorWriter(opt.tumorOpt, (! opt.chromDepthFilename.empty()),
                opt.referenceFilename, cset, tumfs.getStream(),
                opt.isOutputContig),
    rnaWriter(opt.referenceFilename, cset, rnafs.getStream(),
              opt.isOutputContig)
{
    if (0 == opt.edgeOpt.binIndex)
    {
        std::vector<std::string> noSampleNames;
        candWriter.writeHeader(progName, progVersion,noSampleNames);

        const std::vector<std::string>& sampleNames(svScore.sampleNames());
        if (isTumorOnly)
        {
            tumorWriter.writeHeader(progName, progVersion,sampleNames);
        }
        else if (opt.isRNA)
        {
            rnaWriter.writeHeader(progName, progVersion, sampleNames);
        }
        else
        {
            std::vector<std::string> diploidSampleNames(sampleNames.begin(),sampleNames.begin()+svScore.diploidSampleCount());
            diploidWriter.writeHeader(progName, progVersion,diploidSampleNames);
            if (isSomatic) somWriter.writeHeader(progName, progVersion,sampleNames);
        }
    }

    //get_bam_header_sample_name
}



static
bool
isAnyFalse(
    const std::vector<bool>& vb)
{
    for (const bool val : vb)
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
    const std::vector<bool>& isInputJunctionFiltered,
    SupportSamples& svSupports)
{
    const unsigned junctionCount(mjSV.junction.size());

    // relax min spanning count to just 2 for the special case of a junction that is part of a multi-junction EVENT
    //
    const unsigned minJunctionCandidateSpanningCount(std::min(2u,opt.minCandidateSpanningCount));

    // track filtration for each junction:
    std::vector<bool> isJunctionFiltered(isInputJunctionFiltered);

    // first step of SV filtering (before candidates are written):
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
        log_os << __FUNCTION__ << ": Pairs:" << sv.bp1.getPairCount() << " Spanning:" << sv.bp1.getSpanningCount() << "\n";
#endif

        // junction dependent tests:
        //   (1) at least one junction in the set must have spanning count of minCandidateSpanningCount or more
        bool isJunctionSpanFail(false);
        if (isCandidateSpanning)
        {
            if (sv.getPostAssemblySpanningCount(opt.isRNA) < opt.minCandidateSpanningCount)
            {
                isJunctionSpanFail=true;
            }
        }
        if (! isJunctionSpanFail) isCandidateSpanFail=false;

        // independent tests -- as soon as one of these fails, we can continue:
        //   (1) each spanning junction in the set must have spanning count of
        //      minJunctionCandidateSpanningCount or more
        //   (2) no unassembled non-spanning candidates
        if (isCandidateSpanning)
        {
            if (sv.getPostAssemblySpanningCount(opt.isRNA) < minJunctionCandidateSpanningCount)
            {
                isJunctionFiltered[junctionIndex] = true;
                continue;
            }
        }
        else
        {
            if (sv.isImprecise())
            {
                // in this case a non-spanning low-res candidate went into assembly but
                // did not produce a successful contig alignment:
#ifdef DEBUG_GSV
                log_os << __FUNCTION__ << ": Rejecting candidate junction: imprecise non-spanning SV\n";
#endif
                isJunctionFiltered[junctionIndex] = true;
                continue;
            }
        }

        // check min size, or if it is IMPRECISE case where CIEND is a subset of CIPOS for candidate output:
        if (isSVBelowMinSize(sv,opt.scanOpt.minCandidateVariantSize))
        {
#ifdef DEBUG_GSV
            log_os << __FUNCTION__ << ": Filtering out candidate below min size before candidate output stage\n";
#endif
            isJunctionFiltered[junctionIndex] = true;
            continue;
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

        _idgen.getId(edge, sv, opt.isRNA, svId);

        candWriter.writeSV(svData, assemblyData, sv, svId);
    }

    if (opt.isSkipScoring) return;

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

    unsigned sampleCount(svScore.sampleCount());
    unsigned diploidSampleCount(svScore.diploidSampleCount());
    SVModelScoreInfo mjJointModelScoreInfo;
    mjJointModelScoreInfo.setSampleCount(sampleCount, diploidSampleCount);
    svScore.scoreSV(svData, mjAssemblyData, mjSV, junctionSVId,
                    isJunctionFiltered, isSomatic, isTumorOnly,
                    mjModelScoreInfo, mjJointModelScoreInfo,
                    isMJEvent, svSupports);

    const unsigned unfilteredJunctionCount(std::count(isJunctionFiltered.begin(),isJunctionFiltered.end(),true));

    // setup all event-level info that we need to share across all event junctions:
    bool isMJDiploidEvent(isMJEvent);
    EventInfo event;
    event.junctionCount=unfilteredJunctionCount;

    bool isMJEventWriteDiploid(false);
    bool isMJEventWriteSomatic(false);

    std::vector<bool> isJunctionSampleCheckFail(diploidSampleCount, false);

    if (isMJEvent)
    {
        // sample specific check for the assumption of multi-junction candidates: every junction shares the same genotype.
        // check if the individual junctions should be potentially assigned different genotypes.
        for (unsigned sampleIndex(0); sampleIndex<diploidSampleCount; ++sampleIndex)
        {
            const SVScoreInfoDiploid& jointDiploidInfo(mjJointModelScoreInfo.diploid);
            const SVScoreInfoDiploidSample& jointDiploidSampleInfo(jointDiploidInfo.samples[sampleIndex]);
            const DIPLOID_GT::index_t jointGT(jointDiploidSampleInfo.gt);
            const double jointGtPProb(jointDiploidSampleInfo.pprob[jointGT]);

            if (jointGT == DIPLOID_GT::REF)
            {
                isJunctionSampleCheckFail[sampleIndex] = true;
#ifdef DEBUG_GSV
                if (isJunctionSampleCheckFail[sampleIndex])
                {
                    log_os << __FUNCTION__ << ": The multi-junction candidate has hom-ref for sample #"
                           << sampleIndex << ".\n";
                }
#endif
                continue;
            }

            for (unsigned junctionIndex(0); junctionIndex < junctionCount; ++junctionIndex)
            {
                if (isJunctionFiltered[junctionIndex]) continue;

                const SVScoreInfoDiploid& diploidInfo(mjModelScoreInfo[junctionIndex].diploid);
                const SVScoreInfoDiploidSample& diploidSampleInfo(diploidInfo.samples[sampleIndex]);
                const DIPLOID_GT::index_t singleGT(diploidSampleInfo.gt);
                const double singleGtPProb(diploidSampleInfo.pprob[singleGT]);
                const double deltaGtPProb(jointGtPProb - diploidSampleInfo.pprob[jointGT]);

                // fail the genotype check if
                // 1) the joint event changes the genotype of the junction and increases the posterior prob by more than 0.9
                // 2) the posterior prob is larger than 0.9 when the junction is genotyped individually
                if ((! (jointGT==singleGT)) && (deltaGtPProb > 0.9) && (singleGtPProb > 0.9))
                {
                    isJunctionSampleCheckFail[sampleIndex] = true;
                    break;
                }
            }

#ifdef DEBUG_GSV
            if (isJunctionSampleCheckFail[sampleIndex])
            {
                log_os << __FUNCTION__ << ": The multi-junction candidate failed the check of sample-specific genotyping for sample #"
                       << sampleIndex << ".\n";
            }
#endif
        }

        // if any sample passing the genotype check, the multi-junction diploid event remains valid;
        // otherwise, fail the multi-junction candidate
        if (! isAnyFalse(isJunctionSampleCheckFail)) isMJDiploidEvent = false;

#ifdef DEBUG_GSV
        if (! isMJDiploidEvent)
        {
            log_os << __FUNCTION__ << ": Rejecting the multi-junction candidate, failed the sample-specific check of genotyping for all diploid samples.\n";
        }
#endif

        for (unsigned junctionIndex(0); junctionIndex<junctionCount; ++junctionIndex)
        {
            if (isJunctionFiltered[junctionIndex]) continue;

            // set the Id of the first junction in the group as the event (i.e. multi-junction) label
            if (event.label.empty())
            {
                const SVId& svId(junctionSVId[junctionIndex]);
                event.label = svId.localId;
            }

            const SVModelScoreInfo& modelScoreInfo(mjModelScoreInfo[junctionIndex]);

            // for diploid case only,
            // we decide to use multi-junction or single junction based on best score:
            // (for somatic case a lower somatic score could be due to reference evidence in an event member)
            //
            if (mjJointModelScoreInfo.diploid.filters.size() > modelScoreInfo.diploid.filters.size())
            {
                isMJDiploidEvent=false;
            }
            else if (mjJointModelScoreInfo.diploid.altScore < modelScoreInfo.diploid.altScore)
            {
                isMJDiploidEvent=false;
            }

            // for somatic case,
            // report multi-junction if either the multi-junction OR ANY single junction has a good score.
            if ((mjJointModelScoreInfo.somatic.somaticScore >= opt.somaticOpt.minOutputSomaticScore) ||
                (modelScoreInfo.somatic.somaticScore >= opt.somaticOpt.minOutputSomaticScore))
            {
                isMJEventWriteSomatic = true;
            }

            //TODO: set up criteria for isMJEventWriteTumor
        }

        // for events, we write all junctions, or no junctions,
        // so we need to determine write status over the whole set rather than a single junction
        if (isMJDiploidEvent)
        {
            isMJEventWriteDiploid = (mjJointModelScoreInfo.diploid.altScore >= opt.diploidOpt.minOutputAltScore);
        }
    }

    // final scored output is treated (mostly) independently for each junction:
    //
    for (unsigned junctionIndex(0); junctionIndex<junctionCount; ++junctionIndex)
    {
        if (isJunctionFiltered[junctionIndex]) continue;

        const SVCandidateAssemblyData& assemblyData(mjAssemblyData[junctionIndex]);
        const SVCandidate& sv(mjSV.junction[junctionIndex]);
        const SVModelScoreInfo& modelScoreInfo(mjModelScoreInfo[junctionIndex]);

        const SVId& svId(junctionSVId[junctionIndex]);
        const SVScoreInfo& baseInfo(modelScoreInfo.base);
        static const EventInfo nonEvent;

        if (isTumorOnly)
        {
            //TODO: add logic for MJEvent

            const SVScoreInfoTumor& tumorInfo(modelScoreInfo.tumor);
            tumorWriter.writeSV(svData, assemblyData, sv, svId, baseInfo, tumorInfo, nonEvent);
        }
        else if (opt.isRNA)
        {
            const SVScoreInfoRna& rnaInfo(modelScoreInfo.rna);
            rnaWriter.writeSV(svData, assemblyData, sv, svId, baseInfo, rnaInfo, nonEvent);
        }
        else
        {
            {
                const EventInfo& diploidEvent( isMJDiploidEvent ? event : nonEvent );
                const SVModelScoreInfo& scoreInfo( isMJDiploidEvent ? mjJointModelScoreInfo : modelScoreInfo);
                SVScoreInfoDiploid diploidInfo(scoreInfo.diploid);

                if (isMJDiploidEvent)
                {
                    // if one sample failed the genotype check,
                    // use the sample-specific diploid scoreInfo, instead of the joint scoreInfo
                    for (unsigned sampleIndex(0); sampleIndex<diploidSampleCount; ++sampleIndex)
                    {
                        if (isJunctionSampleCheckFail[sampleIndex])
                        {
#ifdef DEBUG_GSV
                            log_os << __FUNCTION__ << ": Junction #" << junctionIndex
                                   << ": Swapped diploid info for sample #" << sampleIndex << ".\n"
                                   << "Before:" << diploidInfo.samples[sampleIndex] << "\n"
                                   << "After:" << modelScoreInfo.diploid.samples[sampleIndex] << "\n" ;
#endif
                            diploidInfo.samples[sampleIndex] = modelScoreInfo.diploid.samples[sampleIndex];
                        }
                    }
                }

                bool isWriteDiploid(false);
                if (isMJDiploidEvent)
                {
                    isWriteDiploid = isMJEventWriteDiploid;
                }
                else
                {
                    isWriteDiploid = (modelScoreInfo.diploid.altScore >= opt.diploidOpt.minOutputAltScore);
                }

                if (isWriteDiploid)
                {
                    diploidWriter.writeSV(svData, assemblyData, sv, svId, baseInfo, diploidInfo, diploidEvent, modelScoreInfo.diploid);
                }
            }

            if (isSomatic)
            {
                const EventInfo& somaticEvent( isMJEvent ? event : nonEvent );
                const SVModelScoreInfo& scoreInfo(isMJEvent ? mjJointModelScoreInfo : modelScoreInfo);
                const SVScoreInfoSomatic& somaticInfo(scoreInfo.somatic);

                bool isWriteSomatic(false);

                if (isMJEvent)
                {
                    isWriteSomatic = isMJEventWriteSomatic;
                }
                else
                {
                    isWriteSomatic = (modelScoreInfo.somatic.somaticScore >= opt.somaticOpt.minOutputSomaticScore);
                }

                if (isWriteSomatic)
                {
                    somWriter.writeSV(svData, assemblyData, sv, svId, baseInfo, somaticInfo, somaticEvent, modelScoreInfo.somatic);
                }
            }
        }
    }
}



SVCandidateProcessor::
SVCandidateProcessor(
    const GSCOptions& opt,
    const SVLocusScanner& readScanner,
    const char* progName,
    const char* progVersion,
    const SVLocusSet& cset,
    EdgeRuntimeTracker& edgeTracker,
    GSCEdgeStatsManager& edgeStatMan) :
    _opt(opt),
    _cset(cset),
    _edgeTracker(edgeTracker),
    _edgeStatMan(edgeStatMan),
    _svRefine(opt, cset.header, cset.getCounts(), _edgeTracker),
    _svWriter(opt, readScanner, cset, progName, progVersion)
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
    _edgeTracker.addCand(isComplex);

    _edgeStatMan.updateJunctionCandidateCounts(edge, junctionCount, isComplex);

    // true if any junction returns a complex (indel) assembly result
    //
    // if true, and if the candidate is multi-junction then score and write each junction independently:
    bool isAnySmallAssembler(false);
    std::vector<SVCandidateAssemblyData> mjAssemblyData(junctionCount);

    if (! _opt.isSkipAssembly)
    {
        const TimeScoper assmTime(_edgeTracker.assmTime);
        for (unsigned junctionIndex(0); junctionIndex<junctionCount; ++junctionIndex)
        {
            const SVCandidate& candidateSV(mjCandidateSV.junction[junctionIndex]);
            SVCandidateAssemblyData& assemblyData(mjAssemblyData[junctionIndex]);
            try
            {
                _svRefine.getCandidateAssemblyData(candidateSV, svData, _opt.isRNA, isFindLargeInsertions, assemblyData);
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
                _edgeTracker.addAssm(isComplex);
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
        const TimeScoper scoreTime(_edgeTracker.scoreTime);
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
