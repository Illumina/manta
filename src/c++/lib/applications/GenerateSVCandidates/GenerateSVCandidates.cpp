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
#include "blt_util/samtools_fasta_util.hh"
#include "common/Exceptions.hh"
#include "common/OutStream.hh"
#include "manta/ReadGroupStatsSet.hh"
#include "manta/SVLocusAssembler.hh"
#include "format/VcfWriterCandidateSV.hh"
#include "format/VcfWriterSomaticSV.hh"
#include "alignment/GlobalJumpAligner.hh"

#include "boost/foreach.hpp"

#include <iostream>
#include <memory>



static
void
getIntervalReferenceSegment(
    const std::string& referenceFilename,
    const bam_header_info& header,
    const pos_t extraRefEdgeSize,
    const GenomeInterval& interval,
    reference_contig_segment& intervalRef)
{
    const bam_header_info::chrom_info& chromInfo(header.chrom_data[interval.tid]);
    const std::string& chrom(chromInfo.label);
    const pos_t beginPos(std::max(0, (interval.range.begin_pos()-extraRefEdgeSize)));
    const pos_t endPos(std::min(static_cast<pos_t>(chromInfo.length), (interval.range.end_pos()+extraRefEdgeSize)));

    // get REF
    intervalRef.set_offset(beginPos);
    get_standardized_region_seq(referenceFilename,chrom,beginPos,endPos,intervalRef.seq());
}



/// extract the reference sequence around each breakend into a reference_contig_segment
/// object
///
/// for each region, we extract the hypothetical breakend region + extraRefEdgeSize bases
/// on each side
///
static
void
getSVReferenceSegments(
    const std::string& referenceFilename,
    const bam_header_info& header,
    const pos_t extraRefEdgeSize,
    const SVCandidate& sv,
    reference_contig_segment& bp1ref,
    reference_contig_segment& bp2ref)
{
    getIntervalReferenceSegment(referenceFilename,header,extraRefEdgeSize,sv.bp1.interval,bp1ref);
    getIntervalReferenceSegment(referenceFilename,header,extraRefEdgeSize,sv.bp2.interval,bp2ref);
}



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



/// we can either traverse all edges in a single locus (disjount subgraph) of the graph
/// OR
/// traverse all edges in out "bin" -- that is, 1 of binCount subsets of the total graph edges
/// designed to be of roughly equal size for parallel processing
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
    int jumpScore(3);
    GlobalJumpAligner<int> aligner(AlignmentScores<int>(1,2,6,0,3),jumpScore);
    while (edger.next())
    {
        const EdgeInfo& edge(edger.getEdge());

        try
        {
            // find number, type and breakend range of SVs on this edge:
            svFind.findCandidateSV(edge,svData,svs);
            BOOST_FOREACH(const SVCandidate& sv, svs)
            {
                Assembly as;
                svAssembler.assembleSVBreakends(sv.bp1,sv.bp2,as);
                svData.getAssembly().insert(svData.getAssembly().end(),
                                            as.begin(),
                                            as.end()
                                            );

                // how much additional reference sequence should we extract from around
                // each side of the breakend region?
                static const pos_t extraRefEdgeSize(200);

                reference_contig_segment bp1ref,bp2ref;
                const std::string bp1RefStr(bp1ref.seq());
                const std::string bp2RefStr(bp2ref.seq());
                getSVReferenceSegments(opt.referenceFilename, cset.header, extraRefEdgeSize, sv, bp1ref, bp2ref);
                for(Assembly::const_iterator ct = as.begin();
                	ct != as.end();
                	++ct) {
                	JumpAlignmentResult<int> res;
                	aligner.align(ct->seq.begin(),ct->seq.end(),
                				  bp1RefStr.begin(),bp1RefStr.end(),
                				  bp2RefStr.begin(),bp2RefStr.end(),
                				  res);
                	svData.getAlignments().push_back(res);
                }
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
