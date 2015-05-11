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
/// \author Xiaoyu Chen
///

#include "format/VcfWriterTumorSV.hh"

void
VcfWriterTumorSV::
addHeaderFormatSampleKey() const
{
    // TODO: extract sample name from input bam header / user
    _os << "\tFORMAT\tSAMPLE";
}


void
VcfWriterTumorSV::
addHeaderInfo() const
{
	_os << "##INFO=<ID=BND_DEPTH,Number=1,Type=Integer,Description=\"Read depth at local translocation breakend\">\n";
	_os << "##INFO=<ID=MATE_BND_DEPTH,Number=1,Type=Integer,Description=\"Read depth at remote translocation mate breakend\">\n";
}

void
VcfWriterTumorSV::
addHeaderFormat() const
{
	_os << "##FORMAT=<ID=PR,Number=.,Type=Integer,Description=\"Spanning paired-read support for the ref and alt alleles in the order listed\">\n";
	_os << "##FORMAT=<ID=SR,Number=.,Type=Integer,Description=\"Split reads for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999\">\n";
}

void
VcfWriterTumorSV::
addHeaderFilters() const
{
	if (_isMaxDepthFilter)
	{
	   _os << "##FILTER=<ID=" << _tumorOpt.maxDepthFilterLabel << ",Description=\"Sample site depth is greater than " << _tumorOpt.maxDepthFactor << "x the mean chromosome depth near one or both variant breakends\">\n";
	}

	_os << "##FILTER=<ID=" << _tumorOpt.maxMQ0FracLabel << ",Description=\"For a small variant (<1000 base), the fraction of reads with MAPQ0 around either breakend exceeds " << _tumorOpt.maxMQ0Frac << "\">\n";
}

void
VcfWriterTumorSV::
writeFilter() const
{
    writeFilters(getTumorInfo().filters);
}


void
VcfWriterTumorSV::
modifySample(
    const SVCandidate& sv,
    SampleTag_t& sampletags) const
{
    const SVScoreInfo& baseInfo(getBaseInfo());

    std::vector<std::string> values(1);

    static const std::string pairTag("PR");
    values[0] = str( boost::format("%i,%i") % baseInfo.normal.ref.confidentSpanningPairCount % baseInfo.normal.alt.confidentSpanningPairCount);
    sampletags.push_back(std::make_pair(pairTag,values));

    if (sv.isImprecise()) return;

    static const std::string srTag("SR");
    values[0] = str( boost::format("%i,%i") % baseInfo.normal.ref.confidentSplitReadCount % baseInfo.normal.alt.confidentSplitReadCount);
    sampletags.push_back(std::make_pair(srTag,values));
}

void
VcfWriterTumorSV::
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
VcfWriterTumorSV::
writeSV(
    const SVCandidateSetData& svData,
    const SVCandidateAssemblyData& adata,
    const SVCandidate& sv,
    const SVId& svId,
    const SVScoreInfo& baseInfo,
    const SVScoreInfoTumor& tumorInfo,
    const EventInfo& event
    )
{
    //TODO: this is a lame way to customize subclass behavior:
    setScoreInfo(baseInfo);
    _tumorInfoPtr=&tumorInfo;
    writeSVCore(svData, adata, sv, svId, event);

    clearScoreInfo();
    _tumorInfoPtr=nullptr;
}
