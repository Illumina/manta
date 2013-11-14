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
    _os << "##FILTER=<ID=" << _diploidOpt.minGTFilterLabel << ",Description=\"GQ score is less than " << _diploidOpt.minGTScoreFilter << "\">\n";
}



#if 0
void
VcfWriterDiploidSV::
modifyInfo(
    InfoTag_t& infotags) const
{
    assert(_modelScorePtr != NULL);
    const SVModelScoreInfo& modelScoreInfo(*_modelScorePtr);

//    infotags.push_back("SOMATIC");
//    infotags.push_back( str(boost::format("SOMATICSCORE=%i") % modelScoreInfo.somatic.somaticScore) );
}
#endif



void
VcfWriterDiploidSV::
modifyTranslocInfo(
    const bool isFirstOfPair,
    InfoTag_t& infotags) const
{
    assert(_modelScorePtr != NULL);
    const SVScoreInfo& baseInfo(_modelScorePtr->base);

    infotags.push_back( str(boost::format("BND_DEPTH=%i") %
                            (isFirstOfPair ? baseInfo.bp1MaxDepth : baseInfo.bp2MaxDepth) ) );
    infotags.push_back( str(boost::format("MATE_BND_DEPTH=%i") %
                            (isFirstOfPair ? baseInfo.bp2MaxDepth : baseInfo.bp1MaxDepth) ) );
}



void
VcfWriterDiploidSV::
writeQual() const
{
    assert(_modelScorePtr != NULL);
    const SVScoreInfoDiploid& diploidInfo(_modelScorePtr->diploid);

    _os << diploidInfo.altScore;
}



void
VcfWriterDiploidSV::
writeFilter() const
{
    assert(_modelScorePtr != NULL);
    const SVScoreInfoDiploid& diploidInfo(_modelScorePtr->diploid);

    writeFilters(diploidInfo.filters);
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
    assert(_modelScorePtr != NULL);
    const SVModelScoreInfo& modelScoreInfo(*_modelScorePtr);
    const SVScoreInfo& baseInfo(modelScoreInfo.base);

    std::vector<std::string> values(1);

    static const std::string gtTag("GT");
    values[0] = gtLabel(modelScoreInfo.diploid.gt);
    sampletags.push_back(std::make_pair(gtTag,values));

    static const std::string gqTag("GQ");
    values[0] = str( boost::format("%s") % modelScoreInfo.diploid.gtScore);
    sampletags.push_back(std::make_pair(gqTag,values));

    static const std::string pairTag("PR");
    values[0] = str( boost::format("%i,%i") % baseInfo.normal.ref.confidentSpanningPairCount % baseInfo.normal.alt.confidentSpanningPairCount);
    sampletags.push_back(std::make_pair(pairTag,values));

    if (sv.isImprecise()) return;

    static const std::string srTag("SR");
    values[0] = str( boost::format("%i,%i") % baseInfo.normal.ref.confidentSplitReadCount % baseInfo.normal.alt.confidentSplitReadCount);
    sampletags.push_back(std::make_pair(srTag,values));
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
