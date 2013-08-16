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

#include "format/VcfWriterSomaticSV.hh"

#include "boost/algorithm/string/join.hpp"

#include "manta/SomaticSVScoreInfo.hh"


void
VcfWriterSomaticSV::
addHeaderInfo() const
{
    _os << "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">\n";
    _os << "##INFO=<ID=SOMATICSCORE,Number=1,Type=Integer,Description=\"Somatic variant Quality score\">\n";
}



void
VcfWriterSomaticSV::
addHeaderFilters() const
{
    if (_isMaxDepthFilter)
    {
        _os << "##FILTER=<ID=" << _somaticOpt.maxDepthFilterLabel << ",Description=\"Normal sample site depth is greater than " << _somaticOpt.maxDepthFactor << "x the mean chromosome depth near one or both variant breakends\">\n";
    }
}



void
VcfWriterSomaticSV::
modifyInfo(
    const SVBreakend& bp1,
    const SVBreakend& bp2,
    const bool isFirstOfPair,
    const SVCandidateSetData& /*svData*/,
    const SVCandidateAssemblyData& adata,
    std::vector<std::string>& infotags) const
{
    assert(_ssInfoPtr != NULL);
    const SomaticSVScoreInfo& ssInfo(*_ssInfoPtr);

    infotags.push_back("SOMATIC");
    infotags.push_back( str(boost::format("SOMATICSCORE=%i") % ssInfo.somaticScore) );
    infotags.push_back( str(boost::format("NORMAL_PAIR_SUPPORT=%i") % ssInfo.normal.spanPairs) );
    infotags.push_back( str(boost::format("TUMOR_PAIR_SUPPORT=%i") % ssInfo.tumor.spanPairs) );
    infotags.push_back( str(boost::format("NORMAL_BND_PAIR_SUPPORT=%i") %
                            (isFirstOfPair ? ssInfo.normal.bp1SpanReads : ssInfo.normal.bp2SpanReads) ) );
    infotags.push_back( str(boost::format("TUMOR_BND_PAIR_SUPPORT=%i") %
                            (isFirstOfPair ? ssInfo.tumor.bp1SpanReads : ssInfo.tumor.bp2SpanReads) ) );
    infotags.push_back( str(boost::format("BND_DEPTH=%i") %
                            (isFirstOfPair ? ssInfo.bp1MaxDepth : ssInfo.bp2MaxDepth) ) );
    infotags.push_back( str(boost::format("MATE_BND_DEPTH=%i") %
                            (isFirstOfPair ? ssInfo.bp2MaxDepth : ssInfo.bp1MaxDepth) ) );

    if (isFirstOfPair)
    {
        // store alignment start + cigar string for each section of the jumping alignment.
        // there can be several contigs per breakend, so we iterate over all of them.
        // only the first breakpoint gets the alignments attached to its VCF entry
        const unsigned numAlign(adata.alignments.size());
        std::string cigar1;
        std::string cigar2;
        for (unsigned alignIndex(0); alignIndex<numAlign; ++alignIndex)
        {
            const SVCandidateAssemblyData::JumpAlignmentResultType align(adata.alignments[alignIndex]);
            infotags.push_back( str(boost::format("CTG_JALIGN_%i_POS_A=%d") %
                                    alignIndex %
                                    (bp1.interval.range.begin_pos()+align.align1.beginPos)) );
            infotags.push_back( str(boost::format("CTG_JALIGN_%i_POS_B=%d") %
                                    alignIndex %
                                    (bp2.interval.range.begin_pos()+align.align2.beginPos)) );

            apath_to_cigar(align.align1.apath,cigar1);
            apath_to_cigar(align.align2.apath,cigar2);

            infotags.push_back( str(boost::format("CTG_JALIGN_%i_CIGAR_A=%s") % alignIndex % cigar1) );
            infotags.push_back( str(boost::format("CTG_JALIGN_%i_CIGAR_B=%s") % alignIndex % cigar2) );
        }
    }
}



std::string
VcfWriterSomaticSV::
getFilter() const
{
    assert(_ssInfoPtr != NULL);
    const SomaticSVScoreInfo& ssInfo(*_ssInfoPtr);

    if (ssInfo.filters.empty())
    {
        return "PASS";
    }
    else
    {
        return boost::algorithm::join(ssInfo.filters, ";");
    }
}


void
VcfWriterSomaticSV::
writeSV(
    const EdgeInfo& edge,
    const SVCandidateSetData& svData,
    const SVCandidateAssemblyData& adata,
    const unsigned svIndex,
    const SVCandidate& sv,
    const SomaticSVScoreInfo& ssInfo)
{

    if (ssInfo.somaticScore < _somaticOpt.minOutputSomaticScore) return;

    //TODO: this is a lame way to customize subclass behavior:
    _ssInfoPtr=&ssInfo;
    writeSVCore(edge, svData, adata, svIndex, sv);
    _ssInfoPtr=NULL;
}
