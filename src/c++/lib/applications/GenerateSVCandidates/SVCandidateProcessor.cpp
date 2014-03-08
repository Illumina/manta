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



void
SVWriter::
writeSV(
    const EdgeInfo& edge,
    const SVCandidateSetData& svData,
    const SVCandidateAssemblyData& assemblyData,
    const SVCandidate& sv)
{
    static const unsigned minCandidateSpanningCount(3);

    const bool isCandidateSpanning(assemblyData.isCandidateSpanning);

#ifdef DEBUG_GSV
    static const std::string logtag("SVWriter::writeSV: ");
    log_os << logtag << "isSpanningSV: " <<  isCandidateSpanning << "\n";
#endif

    if (! isCandidateSpanning)
    {
        if (sv.isImprecise())
        {
            // in this case a non-spanning low-res candidate went into assembly but
            // did not produce a successful contig alignment:
#ifdef DEBUG_GSV
            log_os << logtag << "Rejecting candidate: imprecise non-spanning SV\n";
#endif
            _truthTracker.reportOutcome(SVLog::IMPRECISE_NON_SPANNING);
            return;
        }
    }
    else
    {
        if (sv.bp1.getSpanningCount() < minCandidateSpanningCount)
        {
#ifdef DEBUG_GSV
            log_os << logtag << "Rejecting candidate: minCandidateSpanningCount\n";
#endif
            _truthTracker.reportOutcome(SVLog::LOW_SPANNING_COUNT);
            return;
        }
    }

    // check min size for candidate output:
    if (isSVBelowMinSize(sv,opt.scanOpt.minCandidateVariantSize))
    {
#ifdef DEBUG_GSV
        log_os << logtag << "Filtering out candidate below min size before candidate output stage\n";
#endif
        return;
    }

    candWriter.writeSV(edge, svData, assemblyData, sv);

    // check min size for scoring:
    if (isSVBelowMinSize(sv,opt.minScoredVariantSize))
    {
#ifdef DEBUG_GSV
        log_os << logtag << "Filtering out candidate below min size at scoring stage\n";
#endif
        return;
    }

    svScore.scoreSV(svData, assemblyData, sv, isSomatic, modelScoreInfo);

    if (modelScoreInfo.diploid.altScore >= opt.diploidOpt.minOutputAltScore || opt.isRNA) // todo remove after adding RNA scoring
    {
        diploidWriter.writeSV(edge, svData, assemblyData, sv, modelScoreInfo);
    }

    if (isSomatic)
    {
        if (modelScoreInfo.somatic.somaticScore > opt.somaticOpt.minOutputSomaticScore)
        {
            somWriter.writeSV(edge, svData, assemblyData, sv, modelScoreInfo);
            _truthTracker.reportOutcome(SVLog::WRITTEN);
        }
        else
        {
            _truthTracker.reportOutcome(SVLog::LOW_SOMATIC_SCORE);
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
    assert(! mjCandidateSV.junctions.empty());
    const SVCandidate& candidateSV(mjCandidateSV.junctions[0]);

    if (_opt.isVerbose)
    {
        log_os << __FUNCTION__ << ": Starting analysis of low-resolution candidate: " << candidateSV.candidateIndex << "\n";
    }
#ifdef DEBUG_GSV
    log_os << __FUNCTION__ << ": CandidateSV: " << candidateSV << "\n";
#endif


    const bool isComplex(isComplexSV(candidateSV));
    _edgeTracker.addCand(isComplex);

    SVCandidateAssemblyData assemblyData;

    if (! _opt.isSkipAssembly)
    {
        _svRefine.getCandidateAssemblyData(candidateSV, svData, assemblyData, _opt.isRNA);

        if (_opt.isVerbose)
        {
            log_os << __FUNCTION__ << ": Candidate assembly complete. Assembled candidate count: " << assemblyData.svs.size() << "\n";
        }
    }

    if (assemblyData.svs.empty())
    {
#ifdef DEBUG_GSV
        log_os << __FUNCTION__ << ": score and output low-res candidate\n";
#endif
        _svWriter.writeSV(edge, svData, assemblyData, candidateSV);

    }
    else
    {
        _edgeTracker.addAssm(isComplex);

        _truthTracker.reportNumAssembled(assemblyData.svs.size());

        BOOST_FOREACH(const SVCandidate& assembledSV, assemblyData.svs)
        {
            _truthTracker.addAssembledSV();
#ifdef DEBUG_GSV
            log_os << __FUNCTION__ << ": score and output assembly candidate: " << assembledSV << "\n";
#endif
            _svWriter.writeSV(edge, svData, assemblyData, assembledSV);
        }
    }
}
