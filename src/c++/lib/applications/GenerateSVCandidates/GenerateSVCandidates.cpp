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
#include "format/VcfWriterCandidateSV.hh"
#include "format/VcfWriterSomaticSV.hh"

#include "boost/foreach.hpp"

#include <iostream>
#include <memory>



static
void
dumpEdgeInfo(
    const EdgeInfo& edge,
    const SVLocusSet& set,
    std::ostream& os)
{
    os << "Exception caught while processing graph component: " << edge;
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
    const GSCOptions& opt)
{
    if (opt.isLocusIndex)
    {
        return new EdgeRetrieverLocus(set,opt.locusIndex);
    }
    else
    {
        return new EdgeRetrieverBin(set,opt.binCount,opt.binIndex);
    }
}



static
void
runGSC(
    const GSCOptions& opt,
    const char* progName,
    const char* progVersion)
{
    const bool isSomatic(! opt.somaticOutputFilename.empty());

    SVFinder svFind(opt);
    const SVLocusSet& cset(svFind.getSet());

    // maybe these can be an SV refiner option struct?

    // match, mismatch, open ratios taken from bwa:
    static const AlignmentScores<int> spanningAlignScores(2,-8,-12,-1,-1);
    static const int jumpScore(-25);
    SmallAssemblerOptions spanningAssembleOpt;
    spanningAssembleOpt.minContigLength=75;

    SVCandidateAssemblyRefiner svRefine(opt, cset.header, spanningAssembleOpt, spanningAlignScores, jumpScore);

    SVScorer svScore(opt, cset.header);

    std::auto_ptr<EdgeRetriever> edgerPtr(edgeRFactory(cset, opt));
    EdgeRetriever& edger(*edgerPtr);

    OutStream candfs(opt.candidateOutputFilename);
    OutStream somfs(opt.somaticOutputFilename);

    VcfWriterCandidateSV candWriter(opt.referenceFilename,cset,candfs.getStream());
    VcfWriterSomaticSV somWriter(opt.somaticOpt, (! opt.chromDepthFilename.empty()),
                                 opt.referenceFilename,cset,somfs.getStream());

    if (0 == opt.binIndex)
    {
        candWriter.writeHeader(progName, progVersion);
        if (isSomatic) somWriter.writeHeader(progName, progVersion);
    }

    SVCandidateSetData svData;
    std::vector<SVCandidate> svs;
    SomaticSVScoreInfo ssInfo;

    while (edger.next())
    {
        const EdgeInfo& edge(edger.getEdge());

        try
        {
            // find number, type and breakend range of SVs on this edge:
            svFind.findCandidateSV(edge, svData, svs);

            for (unsigned svIndex(0); svIndex<svs.size(); ++svIndex)
            {
                const SVCandidate& sv(svs[svIndex]);

                SVCandidateAssemblyData assemblyData;
                svRefine.getCandidateAssemblyData(sv, svData, assemblyData);

                const SVCandidate& submitSV(assemblyData.isBestAlignment ? assemblyData.sv : sv);

                candWriter.writeSV(edge, svData, assemblyData, svIndex, submitSV);

                if (isSomatic)
                {
                    svScore.scoreSomaticSV(svData, svIndex, submitSV, ssInfo);
                    somWriter.writeSV(edge, svData, assemblyData, svIndex, submitSV, ssInfo);
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
            dumpEdgeInfo(edge,cset,log_os);
            throw;
        }
    }
}



void
GenerateSVCandidates::
runInternal(int argc, char* argv[]) const
{

    GSCOptions opt;

    parseGSCOptions(*this,argc,argv,opt);
    runGSC(opt, name(), version());
}
