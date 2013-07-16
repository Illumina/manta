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



void
VcfWriterSomaticSV::
addHeaderInfo() const
{
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
    const bool isFirstOfPair,
    std::vector<std::string>& infotags) const
{
    assert(_ssInfoPtr != NULL);
    const SomaticSVScoreInfo& ssInfo(*_ssInfoPtr);

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
    const SVCandidateData& ,
    const unsigned svIndex,
    const SVCandidate& sv,
    const SomaticSVScoreInfo& ssInfo)
{

    if (ssInfo.somaticScore < _somaticOpt.minOutputSomaticScore) return;

    //TODO: this is a lame way to customize subclass behavior:
    _ssInfoPtr=&ssInfo;
    writeSVCore(edge, svIndex, sv);
    _ssInfoPtr=NULL;
}
