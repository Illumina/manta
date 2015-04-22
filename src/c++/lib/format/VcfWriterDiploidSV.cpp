// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
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
    _os << "##INFO=<ID=BND_DEPTH,Number=1,Type=Integer,Description=\"Read depth at local translocation breakend\">\n";
    _os << "##INFO=<ID=MATE_BND_DEPTH,Number=1,Type=Integer,Description=\"Read depth at remote translocation mate breakend\">\n";
    _os << "##INFO=<ID=JUNCTION_QUAL,Number=1,Type=Integer,Description=\"If the SV junction is part of an EVENT (ie. a multi-adjacency variant), this field provides the QUAL value for the adjacency in question only\">\n";
}



void
VcfWriterDiploidSV::
addHeaderFormat() const
{
    _os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    _os << "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">\n";
    _os << "##FORMAT=<ID=PR,Number=.,Type=Integer,Description=\"Spanning paired-read support for the ref and alt alleles in the order listed\">\n";
    _os << "##FORMAT=<ID=SR,Number=.,Type=Integer,Description=\"Split reads for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999\">\n";
}



void
VcfWriterDiploidSV::
addHeaderFilters() const
{
    if (_isMaxDepthFilter)
    {
        _os << "##FILTER=<ID=" << _diploidOpt.maxDepthFilterLabel << ",Description=\"Sample site depth is greater than " << _diploidOpt.maxDepthFactor << "x the mean chromosome depth near one or both variant breakends\">\n";
    }
    _os << "##FILTER=<ID=" << _diploidOpt.minGTFilterLabel << ",Description=\"GQ score is less than " << _diploidOpt.minPassGTScore << "\">\n";
    _os << "##FILTER=<ID=" << _diploidOpt.maxMQ0FracLabel << ",Description=\"For a small variant (<1000 base), the fraction of reads with MAPQ0 around either breakend exceeds " << _diploidOpt.maxMQ0Frac << "\">\n";
}



void
VcfWriterDiploidSV::
modifyInfo(
    const EventInfo& event,
    InfoTag_t& infotags) const
{
    if (event.isEvent())
    {
        infotags.push_back( str(boost::format("JUNCTION_QUAL=%i") % getSingleJunctionDiploidInfo().altScore) );
    }
}



void
VcfWriterDiploidSV::
modifyTranslocInfo(
    const SVCandidate& /*sv*/,
    const bool isFirstOfPair,
    InfoTag_t& infotags) const
{
    const SVScoreInfo& baseInfo(getBaseInfo());

    infotags.push_back( str(boost::format("BND_DEPTH=%i") %
                            (isFirstOfPair ? baseInfo.bp1MaxDepth : baseInfo.bp2MaxDepth) ) );
    infotags.push_back( str(boost::format("MATE_BND_DEPTH=%i") %
                            (isFirstOfPair ? baseInfo.bp2MaxDepth : baseInfo.bp1MaxDepth) ) );
}



void
VcfWriterDiploidSV::
writeQual() const
{
    _os << getDiploidInfo().altScore;
}



void
VcfWriterDiploidSV::
writeFilter() const
{
    writeFilters(getDiploidInfo().filters);
}



static
const char*
gtLabel(
    const DIPLOID_GT::index_t id)
{
    using namespace DIPLOID_GT;
    switch (id)
    {
    case REF :
        return "0/0";
    case HET :
        return "0/1";
    case HOM :
        return "1/1";
    default :
        return "";
    }
}



void
VcfWriterDiploidSV::
modifySample(
    const SVCandidate& sv,
    SampleTag_t& sampletags) const
{
    const SVScoreInfo& baseInfo(getBaseInfo());
    const SVScoreInfoDiploid& diploidInfo(getDiploidInfo());

    std::vector<std::string> values(1);

    static const std::string gtTag("GT");
    values[0] = gtLabel(diploidInfo.gt);
    sampletags.push_back(std::make_pair(gtTag,values));

    static const std::string gqTag("GQ");
    values[0] = str( boost::format("%s") % diploidInfo.gtScore);
    sampletags.push_back(std::make_pair(gqTag,values));

    static const std::string pairTag("PR");
    values[0] = str( boost::format("%i,%i") % baseInfo.normal.ref.confidentSpanningPairCount % baseInfo.normal.alt.confidentSpanningPairCount);
    sampletags.push_back(std::make_pair(pairTag,values));

    if (sv.isImprecise()) return;

    static const std::string srTag("SR");
    values[0] = str( boost::format("%i,%i") % baseInfo.normal.ref.confidentSplitReadCount % baseInfo.normal.alt.confidentSplitReadCount);
    sampletags.push_back(std::make_pair(srTag,values));
    if (_isRNA)
    {
        static const std::string fsTag("FS");
        values[0] = str( boost::format("%i,%i") % baseInfo.normal.ref.splitReadCount % baseInfo.normal.alt.splitReadCount);
        sampletags.push_back(std::make_pair(fsTag,values));
        static const std::string fpTag("FP");
        values[0] = str( boost::format("%i,%i") % baseInfo.normal.ref.spanningPairCount % baseInfo.normal.alt.spanningPairCount);
        sampletags.push_back(std::make_pair(fpTag,values));
    }
}



void
VcfWriterDiploidSV::
writeSV(
    const SVCandidateSetData& svData,
    const SVCandidateAssemblyData& adata,
    const SVCandidate& sv,
    const SVId& svId,
    const SVScoreInfo& baseInfo,
    const SVScoreInfoDiploid& diploidInfo,
    const EventInfo& event,
    const SVScoreInfoDiploid& singleJunctionDiploidInfo)
{
    //TODO: this is a lame way to customize subclass behavior:
    setScoreInfo(baseInfo);
    _diploidInfoPtr=&diploidInfo;
    _singleJunctionDiploidInfoPtr=&singleJunctionDiploidInfo;

    writeSVCore(svData, adata, sv, svId, event);

    clearScoreInfo();
    _diploidInfoPtr=nullptr;
    _singleJunctionDiploidInfoPtr=nullptr;
}
