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
#include "SVFinder.hh"
#include "SVScorer.hh"

#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "common/OutStream.hh"
#include "manta/ReadGroupStatsSet.hh"
#include "manta/SVLocusAssembler.hh"
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

    SVLocusAssembler svAssembler(opt);

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

    SVCandidateData svData;
    std::vector<SVCandidate> svs;
    SomaticSVScoreInfo ssInfo;
    while (edger.next())
    {
        const EdgeInfo& edge(edger.getEdge());

        try
        {
            // find number, type and breakend range of SVs on this edge:
            svFind.findCandidateSV(edge,svData,svs);
            BOOST_FOREACH(const SVCandidate& sv, svs)
            {
            	// TODO: what to do with contigs, how to align them etc?
                // assemble across breakpoint 1
                Assembly a1;
                svAssembler.assembleSVBreakend(sv.bp1,a1);
                svData.getAssemblyBreakp1().insert(svData.getAssemblyBreakp1().end(),
                                                   a1.begin(),
                                                   a1.end()
                                                  );

                // assemble across breakpoint 2
                Assembly a2;
                svAssembler.assembleSVBreakend(sv.bp2,a2);
                svData.getAssemblyBreakp2().insert(svData.getAssemblyBreakp2().end(),
                                                   a2.begin(),
                                                   a2.end());


                // TODO: how to get the reference segment against which to align the contig?
            }

            candWriter.writeSV(edge, svData, svs);

            if (isSomatic)
            {
                unsigned svIndex(0);
                BOOST_FOREACH(const SVCandidate& sv, svs)
                {
                    svScore.scoreSomaticSV(svData, svIndex, sv, ssInfo);
                    somWriter.writeSV(edge, svData, svIndex, sv, ssInfo);
                    svIndex++;
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
