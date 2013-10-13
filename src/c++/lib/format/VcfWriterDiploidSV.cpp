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

#include "format/VcfWriterDiploidSV.hh"

#include "boost/algorithm/string/join.hpp"



void
VcfWriterDiploidSV::
addHeaderFormatSampleKey() const
{
    // TODO: extract sample name from input bam header / user
    _os << "\tFORMAT\tSAMPLE";
}



void
VcfWriterDiploidSV::
addHeaderInfo() const
{
//    _os << "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">\n";
}



void
VcfWriterDiploidSV::
addHeaderFilters() const
{
    if (_isMaxDepthFilter)
    {
        _os << "##FILTER=<ID=" << _diploidOpt.maxDepthFilterLabel << ",Description=\"Sample site depth is greater than " << _diploidOpt.maxDepthFactor << "x the mean chromosome depth near one or both variant breakends\">\n";
    }
}



void
VcfWriterDiploidSV::
addSplitReadInfo(
    InfoTag_t& infotags) const
{
    const SVScoreInfo& baseInfo(_modelScorePtr->base);

    infotags.push_back( str(boost::format("ALT_SPLIT_READ=%i") % baseInfo.normal.alt.splitReadCount));
    infotags.push_back( str(boost::format("REF_SPLIT_READ=%i") % baseInfo.normal.ref.splitReadCount));
}



void
VcfWriterDiploidSV::
modifyInfo(
    InfoTag_t& infotags) const
{
    assert(_modelScorePtr != NULL);
    const SVModelScoreInfo& modelScoreInfo(*_modelScorePtr);

//    infotags.push_back("SOMATIC");
//    infotags.push_back( str(boost::format("SOMATICSCORE=%i") % modelScoreInfo.somatic.somaticScore) );

    const SVScoreInfo& baseInfo(modelScoreInfo.base);
    infotags.push_back( str(boost::format("ALT_PAIR_SUPPORT=%i") % baseInfo.normal.alt.spanPairCount) );
    infotags.push_back( str(boost::format("REF_PAIR_SUPPORT=%i") % baseInfo.normal.ref.spanPairCount) );
}



void
VcfWriterDiploidSV::
modifyTranslocInfo(
    const bool isFirstOfPair,
    InfoTag_t& infotags) const
{
    assert(_modelScorePtr != NULL);
    const SVScoreInfo& baseInfo(_modelScorePtr->base);

    infotags.push_back( str(boost::format("ALT_BND_PAIR_SUPPORT=%i") %
                            (isFirstOfPair ? baseInfo.normal.alt.bp1SpanReadCount : baseInfo.normal.alt.bp2SpanReadCount) ) );
    infotags.push_back( str(boost::format("BND_DEPTH=%i") %
                            (isFirstOfPair ? baseInfo.bp1MaxDepth : baseInfo.bp2MaxDepth) ) );
    infotags.push_back( str(boost::format("MATE_BND_DEPTH=%i") %
                            (isFirstOfPair ? baseInfo.bp2MaxDepth : baseInfo.bp1MaxDepth) ) );
}



std::string
VcfWriterDiploidSV::
getFilter() const
{
    assert(_modelScorePtr != NULL);
    const SVScoreInfo& baseInfo(_modelScorePtr->base);

    if (baseInfo.filters.empty())
    {
        return "PASS";
    }
    else
    {
        return boost::algorithm::join(baseInfo.filters, ";");
    }
}



void
VcfWriterDiploidSV::
modifySample(
    SampleTag_t& sampletags) const
{
    assert(_modelScorePtr != NULL);
//    const SVModelScoreInfo& modelScoreInfo(*_modelScorePtr);

    std::vector<std::string> values(1);
    sampletags.push_back(std::make_pair("GT",values));
}


void
VcfWriterDiploidSV::
writeSV(
    const EdgeInfo& edge,
    const SVCandidateSetData& svData,
    const SVCandidateAssemblyData& adata,
    const SVCandidate& sv,
    const SVModelScoreInfo& modelScore)
{
    //TODO: this is a lame way to customize subclass behavior:
    _modelScorePtr=&modelScore;
    writeSVCore(edge, svData, adata, sv);
    _modelScorePtr=NULL;
}
