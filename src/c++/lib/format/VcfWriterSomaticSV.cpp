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

#include "format/VcfWriterSomaticSV.hh"

#include "boost/algorithm/string/join.hpp"



void
VcfWriterSomaticSV::
addHeaderFormatSampleKey() const
{
    /// TODO: get sample name from bam header/user
    _os << "\tFORMAT\tNORMAL\tTUMOR";
}



void
VcfWriterSomaticSV::
addHeaderInfo() const
{
    _os << "##INFO=<ID=BND_DEPTH,Number=1,Type=Integer,Description=\"Read depth at local translocation breakend\">\n";
    _os << "##INFO=<ID=MATE_BND_DEPTH,Number=1,Type=Integer,Description=\"Read depth at remote translocation mate breakend\">\n";
    _os << "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">\n";
    _os << "##INFO=<ID=SOMATICSCORE,Number=1,Type=Integer,Description=\"Somatic variant quality score\">\n";
}



void
VcfWriterSomaticSV::
addHeaderFormat() const
{
    _os << "##FORMAT=<ID=PR,Number=.,Type=Integer,Description=\"Spanning paired-read support for the ref and alt alleles in the order listed\">\n";
    _os << "##FORMAT=<ID=SR,Number=.,Type=Integer,Description=\"Split reads for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999\">\n";
}



void
VcfWriterSomaticSV::
addHeaderFilters() const
{
    if (_isMaxDepthFilter)
    {
        _os << "##FILTER=<ID=" << _somaticOpt.maxDepthFilterLabel << ",Description=\"Normal sample site depth is greater than " << _somaticOpt.maxDepthFactor << "x the mean chromosome depth near one or both variant breakends\">\n";
    }
    _os << "##FILTER=<ID=" << _somaticOpt.minSomaticScoreLabel << ",Description=\"Somatic score is less than " << _somaticOpt.minPassSomaticScore << "\">\n";
    _os << "##FILTER=<ID=" << _somaticOpt.maxMQ0FracLabel << ",Description=\"For a small variant (<1000 bases) in the normal sample, the fraction of reads with MAPQ0 around either breakend exceeds " << _somaticOpt.maxMQ0Frac << "\">\n";
}



void
VcfWriterSomaticSV::
modifyInfo(
    std::vector<std::string>& infotags) const
{
    infotags.push_back("SOMATIC");
    infotags.push_back( str(boost::format("SOMATICSCORE=%i") % getSomaticInfo().somaticScore) );
}



void
VcfWriterSomaticSV::
modifyTranslocInfo(
    const bool isFirstOfPair,
    std::vector<std::string>& infotags) const
{
    const SVScoreInfo& baseInfo(getBaseInfo());

    infotags.push_back( str(boost::format("BND_DEPTH=%i") %
                            (isFirstOfPair ? baseInfo.bp1MaxDepth : baseInfo.bp2MaxDepth) ) );
    infotags.push_back( str(boost::format("MATE_BND_DEPTH=%i") %
                            (isFirstOfPair ? baseInfo.bp2MaxDepth : baseInfo.bp1MaxDepth) ) );
}



void
VcfWriterSomaticSV::
modifySample(
    const SVCandidate& sv,
    SampleTag_t& sampletags) const
{
    const SVScoreInfo& baseInfo(getBaseInfo());

    std::vector<std::string> values(2);

    static const std::string pairTag("PR");
    values[0] = str( boost::format("%i,%i") % baseInfo.normal.ref.confidentSpanningPairCount % baseInfo.normal.alt.confidentSpanningPairCount);
    values[1] = str( boost::format("%i,%i") % baseInfo.tumor.ref.confidentSpanningPairCount % baseInfo.tumor.alt.confidentSpanningPairCount);
    sampletags.push_back(std::make_pair(pairTag,values));

    if (sv.isImprecise()) return;

    static const std::string srTag("SR");
    values[0] = str( boost::format("%i,%i") % baseInfo.normal.ref.confidentSplitReadCount % baseInfo.normal.alt.confidentSplitReadCount);
    values[1] = str( boost::format("%i,%i") % baseInfo.tumor.ref.confidentSplitReadCount % baseInfo.tumor.alt.confidentSplitReadCount);
    sampletags.push_back(std::make_pair(srTag,values));
}



void
VcfWriterSomaticSV::
writeFilter() const
{
    writeFilters(getSomaticInfo().filters);
}



void
VcfWriterSomaticSV::
writeSV(
    const SVCandidateSetData& svData,
    const SVCandidateAssemblyData& adata,
    const SVCandidate& sv,
    const SVId& svId,
    const SVScoreInfo& baseInfo,
    const SVScoreInfoSomatic& somaticInfo,
    const EventInfo& event)
{
    //TODO: this is a lame way to customize subclass behavior:
    setScoreInfo(baseInfo);
    _somaticInfoPtr=&somaticInfo;

    writeSVCore(svData, adata, sv, svId, event);

    clearScoreInfo();
    _somaticInfoPtr=NULL;
}
