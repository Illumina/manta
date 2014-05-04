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

#include "GenerateSVCandidates.hh"
#include "EdgeRetrieverJumpBin.hh"
#include "EdgeRetrieverLocus.hh"
#include "EdgeRuntimeTracker.hh"
#include "GSCOptions.hh"
#include "SVCandidateProcessor.hh"
#include "SVFinder.hh"

#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "manta/MultiJunctionUtil.hh"
#include "manta/SVCandidateUtil.hh"
#include "svgraph/EdgeInfoUtil.hh"
#include "truth/TruthTracker.hh"

#include "boost/foreach.hpp"

#include <iostream>

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
        return (new EdgeRetrieverLocus(set, opt.graphNodeMaxEdgeCount, opt.locusOpt));
    }
    else
    {
        return (new EdgeRetrieverJumpBin(set, opt.graphNodeMaxEdgeCount, opt.binCount, opt.binIndex));
    }
}



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

    TruthTracker truthTracker(opt.truthVcfFilename, cset);
    EdgeRuntimeTracker edgeTracker(opt.edgeRuntimeFilename);

    SVCandidateProcessor svProcessor(opt, progName, progVersion, cset,  truthTracker, edgeTracker);

    std::auto_ptr<EdgeRetriever> edgerPtr(edgeRFactory(cset, opt.edgeOpt));
    EdgeRetriever& edger(*edgerPtr);

    SVCandidateSetData svData;
    std::vector<SVCandidate> svs;
    std::vector<SVMultiJunctionCandidate> mjSVs;

    static const std::string logtag("runGSC");
    if (opt.isVerbose)
    {
        log_os << logtag << " " << cset.header << "\n";
    }

    while (edger.next())
    {
        const EdgeInfo& edge(edger.getEdge());

        try
        {
            truthTracker.addEdge(edge);
            edgeTracker.start();

            if (opt.isVerbose)
            {
                log_os << logtag << " starting analysis of edge: ";
                dumpEdgeInfo(edge,cset,log_os);
            }

            // determine if this is the only edge for this node:
            const bool isIsolatedEdge(testIsolatedEdge(cset,edge));

            // find number, type and breakend range (or better: breakend distro) of SVs on this edge:
            svFind.findCandidateSV(edge, svData, svs, truthTracker);

            truthTracker.reportNumCands(svs.size(), edge);

            if (opt.isVerbose)
            {
                log_os << logtag << " Low-resolution candidate generation complete. Candidate count: " << svs.size() << "\n";
            }

            const unsigned svCount(svs.size());
            for (unsigned i(0); i<svCount; ++i)
            {
                truthTracker.addCandSV();
            }

            findMultiJunctionCandidates(svs, mjSVs);

            bool isFindLargeInsertions(isIsolatedEdge);
            if (isFindLargeInsertions)
            {
                BOOST_FOREACH(const SVMultiJunctionCandidate& mjCandidateSV, mjSVs)
                {
                    BOOST_FOREACH(const SVCandidate& candidateSV, mjCandidateSV.junction)
                    {
                        if (! isComplexSV(candidateSV)) isFindLargeInsertions=false;
                    }
                }
            }

            BOOST_FOREACH(const SVMultiJunctionCandidate& mjCandidateSV, mjSVs)
            {
                svProcessor.evaluateCandidate(edge, mjCandidateSV, svData, isFindLargeInsertions);
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

        edgeTracker.stop(edge);
        if (opt.isVerbose)
        {
            log_os << logtag << " Processing this edge took " << edgeTracker.getLastEdgeTime() << " seconds.\n";
        }
    }

    truthTracker.dumpAll();
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
