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
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "GenerateSVCandidates.hh"
#include "EdgeRetrieverBin.hh"
#include "EdgeRetrieverLocus.hh"
#include "GSCOptions.hh"
#include "SVCandidateAssemblyRefiner.hh"
#include "SVFinder.hh"
#include "SVScorer.hh"

#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "common/OutStream.hh"
#include "manta/ReadGroupStatsSet.hh"
#include "manta/SVCandidateAssemblyData.hh"
#include "manta/SVCandidateUtil.hh"
#include "format/VcfWriterCandidateSV.hh"
#include "format/VcfWriterSomaticSV.hh"

#include "boost/foreach.hpp"

#include <iostream>
#include <memory>
#include <time.h>

//#define DEBUG_GSV


/// provide additional edge details, intended for attachment to an in-flight exception:
static
void
dumpEdgeInfo(
    const EdgeInfo& edge,
    const SVLocusSet& set,
    std::ostream& os)
{
    os << edge;
    os << "\tnode1:" << set.getLocus(edge.locusIndex).getNode(edge.nodeIndex1);
    os << "\tnode2:" << set.getLocus(edge.locusIndex).getNode(edge.nodeIndex2);
}



/// we can either traverse all edges in a single locus (disjoint subgraph) of the graph
/// OR
/// traverse all edges in one "bin" -- that is, one out of binCount subsets of the total
/// graph edges. Each bin is designed to be of roughly equal size in terms of total
/// anticipated workload, so that we have good parallel processing performance.
///
static
EdgeRetriever*
edgeRFactory(
    const SVLocusSet& set,
    const EdgeOptions& opt)
{
    if (opt.isLocusIndex)
    {
        return (new EdgeRetrieverLocus(set, opt.graphNodeMaxEdgeCount, opt.locusIndex));
    }
    else
    {
        return (new EdgeRetrieverBin(set, opt.graphNodeMaxEdgeCount, opt.binCount, opt.binIndex));
    }
}



/// hack object setup to allow iteration over multiple sv candidates
///
struct SVWriter
{
    SVWriter(
        const GSCOptions& initOpt,
        const SVLocusSet& cset,
        const char* progName,
        const char* progVersion) :
        opt(initOpt),
        isSomatic(! opt.somaticOutputFilename.empty()),
        svScore(opt, cset.header),
        candfs(opt.candidateOutputFilename),
        somfs(opt.somaticOutputFilename),
        candWriter(opt.referenceFilename,cset,candfs.getStream()),
        somWriter(opt.somaticOpt, (! opt.chromDepthFilename.empty()),
                  opt.referenceFilename,cset,somfs.getStream())
    {
        if (0 == opt.edgeOpt.binIndex)
        {
            candWriter.writeHeader(progName, progVersion);
            if (isSomatic) somWriter.writeHeader(progName, progVersion);
        }
    }

    void
    writeSV(
        const EdgeInfo& edge,
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& assemblyData,
        const SVCandidate& sv)
    {
        static const unsigned minCandidatePairCount(3);

        const bool isSelfEdge(edge.nodeIndex1 == edge.nodeIndex2);

#ifdef DEBUG_GSV
        static const std::string logtag("SVWriter::writeSV");
        log_os << logtag << " isSelfEdge: " <<  isSelfEdge << "\n";
#endif

        if (isSelfEdge)
        {
            if (sv.isImprecise())
            {
#ifdef DEBUG_GSV
                log_os << logtag << " rejecting candidate, imprecise self-edge\n";
#endif
                return;
            }
        }
        else
        {
            if (sv.bp1.pairCount < minCandidatePairCount)
            {
#ifdef DEBUG_GSV
                log_os << logtag << " rejecting candidate, minCandidatePairCount\n";
#endif

                return;
            }
        }

        // check min size for candidate output:
        if (isSVBelowMinSize(sv,opt.scanOpt.minCandidateVariantSize))
        {
#ifdef DEBUG_GSV
            log_os << logtag << " filtering out candidate below min size before candidate output stage\n";
#endif
            return;
        }

        candWriter.writeSV(edge, svData, assemblyData, sv);

        // check min size for scoring:
        if (isSVBelowMinSize(sv,opt.minScoredVariantSize))
        {
#ifdef DEBUG_GSV
            log_os << logtag << " filtering out candidate below min size at scoring stage\n";
#endif
            return;
        }

        if (isSomatic)
        {
            svScore.scoreSomaticSV(svData, assemblyData, sv, ssInfo);

            if (ssInfo.somaticScore > opt.somaticOpt.minOutputSomaticScore)
            {
                somWriter.writeSV(edge, svData, assemblyData, sv, ssInfo);
            }
        }
    }

    ///////////////////////// data:
    const GSCOptions& opt;
    const bool isSomatic;

    SVScorer svScore;
    SomaticSVScoreInfo ssInfo;

    OutStream candfs;
    OutStream somfs;

    VcfWriterCandidateSV candWriter;
    VcfWriterSomaticSV somWriter;
};


#if 0
/// edge indices+graph evidence counts and regions:
///
/// this is designed to be useful even when the locus graph is not present
struct EhancedEdgeInfo
{

};

/// reduce the full (very-large) graph down to just the information we need during SVCandidate generation:
struct ReducedGraphInfo
{
    ReducedGraphInfo(const GSCOptions& opt)

    bam_header_info header;

    std::vector<EnhancedEdgeInfo> edges;
};
#endif


static
void
runGSC(
    const GSCOptions& opt,
    const char* progName,
    const char* progVersion)
{
#if 0
    {
        // to save memory, load the graph and process/store only the information we need from it:
    }
#endif

    SVFinder svFind(opt);
    const SVLocusSet& cset(svFind.getSet());

    SVCandidateAssemblyRefiner svRefine(opt, cset.header);

    std::auto_ptr<EdgeRetriever> edgerPtr(edgeRFactory(cset, opt.edgeOpt));
    EdgeRetriever& edger(*edgerPtr);

    SVWriter svWriter(opt, cset, progName, progVersion);

    SVCandidateSetData svData;
    std::vector<SVCandidate> svs;

    static const std::string logtag("runGSC");
    if (opt.isVerbose)
    {
        log_os << logtag << " " << cset.header << "\n";
    }

    while (edger.next())
    {
        const EdgeInfo& edge(edger.getEdge());

        clock_t startTime(0);
        if (opt.isVerbose)
        {
            log_os << logtag << " starting analysis of edge: ";
            dumpEdgeInfo(edge,cset,log_os);
            startTime = clock();
        }

        try
        {
            // find number, type and breakend range (or better: breakend distro) of SVs on this edge:
            svFind.findCandidateSV(edge, svData, svs);

            if (opt.isVerbose)
            {
                log_os << logtag << " Low-resolution candidate generation complete. Candidate count: " << svs.size() << "\n";
            }

            BOOST_FOREACH(const SVCandidate& candidateSV, svs)
            {
                if (opt.isVerbose)
                {
                    log_os << logtag << " Starting analysis of low-resolution candidate: " << candidateSV.candidateIndex << "\n";
                }
#ifdef DEBUG_GSV
                log_os << logtag << " CandidateSV: " << candidateSV << "\n";
#endif
                SVCandidateAssemblyData assemblyData;

                if (! opt.isSkipAssembly)
                {
                    svRefine.getCandidateAssemblyData(candidateSV, svData, assemblyData);

                    if (opt.isVerbose)
                    {
                        log_os << logtag << " Candidate assembly complete. Assembled candidate count: " << assemblyData.svs.size() << "\n";
                    }
                }

                if (assemblyData.svs.empty())
                {
#ifdef DEBUG_GSV
                    log_os << logtag << " score and output low-res candidate\n";
#endif
                    svWriter.writeSV(edge, svData, assemblyData, candidateSV);
                }
                else
                {
                    BOOST_FOREACH(const SVCandidate& assembledSV, assemblyData.svs)
                    {
#ifdef DEBUG_GSV
                        log_os << logtag << " score and output assembly candidate: " << assembledSV << "\n";
#endif
                        svWriter.writeSV(edge, svData, assemblyData, assembledSV);
                    }
                }
            }
        }
        catch (illumina::common::ExceptionData& e)
        {
            std::ostringstream oss;
            dumpEdgeInfo(edge,cset,oss);
            e << illumina::common::ExceptionMsg(oss.str());
            throw;
        }
        catch (...)
        {
            log_os << "Exception caught while processing graph component: ";
            dumpEdgeInfo(edge,cset,log_os);
            throw;
        }

        if (opt.isVerbose)
        {
            static const double clockFactor(1./static_cast<double>(CLOCKS_PER_SEC));

            const double elapsedSec(( clock() - startTime )*clockFactor);
            log_os << logtag << " Processing this edge took " << elapsedSec << " seconds.\n";
        }
    }
}



void
GenerateSVCandidates::
runInternal(int argc, char* argv[]) const
{
    GSCOptions opt;

    parseGSCOptions(*this,argc,argv,opt);
#ifdef DEBUG_GSV
    opt.isVerbose=true;
#endif
    runGSC(opt, name(), version());
}
