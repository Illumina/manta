//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
//
//

/// \file
/// \author Xiaoyu Chen
///

#include "format/VcfWriterTumorSV.hpp"

void VcfWriterTumorSV::addHeaderInfo(std::ostream& os) const
{
  os << "##INFO=<ID=BND_DEPTH,Number=1,Type=Integer,Description=\"Read depth at local translocation breakend\">\n";
  os << "##INFO=<ID=MATE_BND_DEPTH,Number=1,Type=Integer,Description=\"Read depth at remote translocation mate breakend\">\n";
}

void VcfWriterTumorSV::addHeaderFormat(std::ostream& os) const
{
  os << "##FORMAT=<ID=PR,Number=.,Type=Integer,Description=\"Spanning paired-read support for the ref and alt alleles in the order listed\">\n";
  os << "##FORMAT=<ID=SR,Number=.,Type=Integer,Description=\"Split reads for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999\">\n";
}

void VcfWriterTumorSV::addHeaderFilters(std::ostream& os) const
{
  if (_isMaxDepthFilter) {
    os << "##FILTER=<ID=" << _tumorOpt.maxDepthFilterLabel
       << ",Description=\"Sample site depth is greater than " << _tumorOpt.maxDepthFactor
       << "x the median chromosome depth near one or both variant breakends\">\n";
  }

  os << "##FILTER=<ID=" << _tumorOpt.maxMQ0FracLabel
     << ",Description=\"For a small variant (<1000 base), the fraction of reads with MAPQ0 around either breakend exceeds "
     << _tumorOpt.maxMQ0Frac << "\">\n";
}

void VcfWriterTumorSV::writeFilter(const boost::any specializedScoringInfo, std::ostream& os) const
{
  const SVScoreInfoTumor& tumorScoringInfo(*boost::any_cast<const SVScoreInfoTumor*>(specializedScoringInfo));
  writeFilters(tumorScoringInfo.filters, os);
}

void VcfWriterTumorSV::modifySample(
    const SVCandidate& sv,
    const SVScoreInfo* baseScoringInfoPtr,
    const boost::any /*specializedScoringInfo*/,
    SampleTag_t& sampletags) const
{
  assert(baseScoringInfoPtr);
  const SVScoreInfo& baseScoringInfo(*baseScoringInfoPtr);

  const unsigned sampleCount(baseScoringInfo.samples.size());

  std::vector<std::string> values(sampleCount);

  for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex) {
    const SVSampleInfo& sinfo(baseScoringInfo.samples[sampleIndex]);
    values[sampleIndex] = str(
        boost::format("%i,%i") % sinfo.ref.confidentSpanningPairCount % sinfo.alt.confidentSpanningPairCount);
  }
  sampletags.push_back(std::make_pair("PR", values));

  if (sv.isImprecise()) return;

  for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex) {
    const SVSampleInfo& sinfo(baseScoringInfo.samples[sampleIndex]);
    values[sampleIndex] =
        str(boost::format("%i,%i") % sinfo.ref.confidentSplitReadCount % sinfo.alt.confidentSplitReadCount);
  }
  sampletags.push_back(std::make_pair("SR", values));
}

void VcfWriterTumorSV::modifyTranslocInfo(
    const SVCandidate& /*sv*/,
    const SVScoreInfo* baseScoringInfoPtr,
    const bool         isFirstOfPair,
    const SVCandidateAssemblyData& /*assemblyData*/,
    InfoTag_t& infotags) const
{
  assert(baseScoringInfoPtr);
  const SVScoreInfo& baseScoringInfo(*baseScoringInfoPtr);
  infotags.push_back(
      str(boost::format("BND_DEPTH=%i") %
          (isFirstOfPair ? baseScoringInfo.bp1MaxDepth : baseScoringInfo.bp2MaxDepth)));
  infotags.push_back(
      str(boost::format("MATE_BND_DEPTH=%i") %
          (isFirstOfPair ? baseScoringInfo.bp2MaxDepth : baseScoringInfo.bp1MaxDepth)));
}

void VcfWriterTumorSV::writeSV(
    const SVCandidateSetData&      svData,
    const SVCandidateAssemblyData& adata,
    const SVCandidate&             sv,
    const SVId&                    svId,
    const SVScoreInfo&             baseScoringInfo,
    const SVScoreInfoTumor&        tumorInfo,
    const EventInfo&               event) const
{
  writeSVCore(svData, adata, sv, svId, &baseScoringInfo, &tumorInfo, event);
}
