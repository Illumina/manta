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
#include "GSCOptions.hh"
#include "EdgeRetriever.hh"
#include "SVFinder.hh"
#include "SVScorer.hh"
#include "WriteCandidateSVToVcf.hh"
#include "WriteSomaticSVToVcf.hh"

#include "common/OutStream.hh"
#include "manta/ReadGroupStatsSet.hh"

#include "boost/foreach.hpp"

#include <iostream>




static
void
runGSC(
    const GSCOptions& opt,
    const char* version)
{
    {
        // early test that we have permission to write to output file
        OutStream outs(opt.candidateOutputFilename);

        if (! opt.somaticOutputFilename.empty())
        {
            OutStream somouts(opt.somaticOutputFilename);
        }
    }

    SVFinder svFind(opt);
    SVScorer svScore(opt);

    // load in set:
    SVLocusSet set;
    set.load(opt.graphFilename.c_str());
    const SVLocusSet& cset(svFind.getSet());

    EdgeRetriever edger(cset, opt.binCount, opt.binIndex);

    OutStream outs(opt.candidateOutputFilename);
    std::ostream& outfp(outs.getStream());

    const bool isSomatic(! opt.somaticOutputFilename.empty());

    if (0 == opt.binIndex)
    {
        writeCandidateSVVcfHeader(opt.referenceFilename.c_str(), version, outfp);

        if(isSomatic)
        {
            writeSomaticSVVcfHeader(opt.referenceFilename.c_str(), version, outfp);
        }
    }

    SVCandidateData svData;
    std::vector<SVCandidate> svs;
    SomaticSVScoreInfo ssInfo;
    while (edger.next())
    {
        const EdgeInfo& edge(edger.getEdge());

        // find number, type and breakend range of SVs on this edge:
        svFind.findCandidateSV(edge,svData,svs);

        writeCandidateSVToVcf(opt.referenceFilename, set, edge, svData, svs, outfp);

        if(isSomatic)
        {
            const unsigned svIndex(0);
            BOOST_FOREACH(const SVCandidate& sv, svs)
            {
                svScore.scoreSomaticSV(svData, svIndex,sv,ssInfo);
                writeSomaticSVToVcf(opt.referenceFilename, set, edge, svData, svIndex, sv, ssInfo, outfp);
            }
        }
    }

}



void
GenerateSVCandidates::
runInternal(int argc, char* argv[]) const
{

    GSCOptions opt;

    parseGSCOptions(*this,argc,argv,opt);
    runGSC(opt,version());
}
