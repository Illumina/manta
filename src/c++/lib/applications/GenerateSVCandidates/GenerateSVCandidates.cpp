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

#include "common/OutStream.hh"
#include "manta/ReadGroupStatsSet.hh"

#include "boost/foreach.hpp"

#include <iostream>



static
void
writeVcfHeader(std::ostream& os)
{
    os << "##fileformat=VCFv4.1\n";
    os << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
}


static
void
writeTransloc(
        const bam_header_info& header,
        const SVBreakend& bp1,
        const SVBreakend& bp2,
        std::ostream& os)
{
    const std::string& chrom1(header.chrom_data[bp1.interval.tid].label);
    const std::string& chrom2(header.chrom_data[bp2.interval.tid].label);

    pos_t pos(0);

    os << chrom1
       << '\t' << pos
       << '\t' << '.' // ID
       << '\t' << '.' // REF
       << '\t' << chrom2 // ALT
       << '\t' << '.' // QUAL
       << '\t' << '.' // FILTER
       << '\t' << '.' // INFO
       << '\n';
}


static
void
writeCandidates(
        const SVLocusSet& set,
        const SVCandidateData& ,
        const std::vector<SVCandidate>& svs,
        std::ostream& outfp)
{
    const unsigned minPairCount(set.getMinMergeEdgeCount());
    const bam_header_info& header(set.header);

    BOOST_FOREACH(const SVCandidate& sv, svs)
    {
        if(sv.bp1.pairCount < minPairCount) continue;

        const SV_TYPE::index_t svType(getSVType(sv));
        if(svType == SV_TYPE::INTERTRANSLOC)
        {
            writeTransloc(header, sv.bp1, sv.bp2, outfp);
            writeTransloc(header, sv.bp2, sv.bp1, outfp);
        }
        else
        {
            assert(! "Unknown SV type");
        }
        //write debug output
        outfp << sv;

        //scoreSV
    }
}


static
void
runGSC(const GSCOptions& opt)
{
    {
        // early test that we have permission to write to output file
        OutStream outs(opt.candidateOutputFilename);

        if (! opt.somaticOutputFilename.empty())
        {
            OutStream somouts(opt.somaticOutputFilename);
        }
    }

    SVFinder finder(opt);

    // load in set:
    SVLocusSet set;
    set.load(opt.graphFilename.c_str());
    const SVLocusSet& cset(finder.getSet());

    EdgeRetriever edger(cset, opt.binCount, opt.binIndex);

    OutStream outs(opt.candidateOutputFilename);
    std::ostream& outfp(outs.getStream());
    writeVcfHeader(outfp);


    SVCandidateData svData;
    std::vector<SVCandidate> svs;
    while (edger.next())
    {
        const EdgeInfo& edge(edger.getEdge());

        // find number, type and breakend range of SVs on this edge:
        finder.findSVCandidates(edge,svData,svs);

        writeCandidates(set,svData,svs,outfp);
    }

}



void
GenerateSVCandidates::
runInternal(int argc, char* argv[]) const
{

    GSCOptions opt;

    parseGSCOptions(*this,argc,argv,opt);
    runGSC(opt);
}
