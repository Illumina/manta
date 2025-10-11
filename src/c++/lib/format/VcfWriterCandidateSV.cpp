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
/// \author Chris Saunders
///

#include "format/VcfWriterCandidateSV.hpp"

void VcfWriterCandidateSV::addHeaderInfo(std::ostream& os) const
{
  os << "##INFO=<ID=PAIR_COUNT,Number=1,Type=Integer,Description=\"Read pairs supporting this variant where both reads are confidently mapped\">\n";
  os << "##INFO=<ID=BND_PAIR_COUNT,Number=1,Type=Integer,Description=\"Confidently mapped reads supporting this variant at this breakend (mapping may not be confident at remote breakend)\">\n";
  os << "##INFO=<ID=UPSTREAM_PAIR_COUNT,Number=1,Type=Integer,Description=\"Confidently mapped reads supporting this variant at the upstream breakend (mapping may not be confident at downstream breakend)\">\n";
  os << "##INFO=<ID=DOWNSTREAM_PAIR_COUNT,Number=1,Type=Integer,Description=\"Confidently mapped reads supporting this variant at this downstream breakend (mapping may not be confident at upstream breakend)\">\n";
}

void VcfWriterCandidateSV::modifyTranslocInfo(
    const SVCandidate& sv,
    const SVScoreInfo* /*baseScoringInfoPtr*/,
    const bool isFirstOfPair,
    const SVCandidateAssemblyData& /*assemblyData*/,
    InfoTag_t& infoTags) const
{
  const SVBreakend& bpA(isFirstOfPair ? sv.bp1 : sv.bp2);

  infoTags.push_back(str(boost::format("BND_PAIR_COUNT=%i") % bpA.getLocalPairCount()));
  infoTags.push_back(str(boost::format("PAIR_COUNT=%i") % bpA.getPairCount()));
}

void VcfWriterCandidateSV::modifyInvdelInfo(
    const SVCandidate& sv, const bool isBp1First, InfoTag_t& infoTags) const
{
  const SVBreakend& bpA(isBp1First ? sv.bp1 : sv.bp2);
  const SVBreakend& bpB(isBp1First ? sv.bp2 : sv.bp1);

  infoTags.push_back(str(boost::format("UPSTREAM_PAIR_COUNT=%i") % bpA.getLocalPairCount()));
  infoTags.push_back(str(boost::format("DOWNSTREAM_PAIR_COUNT=%i") % bpB.getLocalPairCount()));
  infoTags.push_back(str(boost::format("PAIR_COUNT=%i") % bpA.getPairCount()));
}

void VcfWriterCandidateSV::writeSV(
    const SVCandidateSetData&      svData,
    const SVCandidateAssemblyData& adata,
    const SVCandidate&             sv,
    const SVId&                    svId) const
{
  static const EventInfo event;
  writeSVCore(svData, adata, sv, svId, nullptr, nullptr, event);
}
