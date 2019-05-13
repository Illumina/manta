//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Felix Schlesinger
///

#include "format/VcfWriterRnaSV.hpp"

void VcfWriterRnaSV::addHeaderInfo(std::ostream& os) const
{
  os << "##INFO=<ID=BND_DEPTH,Number=1,Type=Integer,Description=\"Read depth at local translocation breakend\">\n";
  os << "##INFO=<ID=MATE_BND_DEPTH,Number=1,Type=Integer,Description=\"Read depth at remote translocation mate breakend\">\n";
  os << "##INFO=<ID=REF_COUNT,Number=1,Type=Integer,Description=\"The number of reads supporting the reference allele at this breakend\">\n";
  os << "##INFO=<ID=MATE_REF_COUNT,Number=1,Type=Integer,Description=\"The number of reads supporting the reference allele at the other breakend\">\n";
  os << "##INFO=<ID=RNA_FIRST,Number=0,Type=Flag,Description=\"For RNA fusions, this break-end is 5' in the fusion transcript\">\n";
  os << "##INFO=<ID=RNA_STRANDED,Number=0,Type=Flag,Description=\"For RNA fusions, the direction of transcription is known\">\n";
  os << "##INFO=<ID=RNA_FwRvReads,Number=2,Type=Integer,Description=\"For RNA fusions, number of stranded reads supporting forward or reverse direction of transcription\">\n";
  os << "##INFO=<ID=RNA_Reads,Number=1,Type=Integer,Description=\"The number of reads and pairs that potentially support this candidate before refinement and scoring\">\n";
  os << "##INFO=<ID=RNA_CONTIG,Number=1,Type=String,Description=\"The sequence of the breakend spanning contig\">\n";
  os << "##INFO=<ID=RNA_CONTIG_ALN,Number=2,Type=Integer,Description=\"Length of the spanning contig alignment on each breakend\">\n";
}

void VcfWriterRnaSV::addHeaderFormat(std::ostream& os) const
{
  os << "##FORMAT=<ID=PR,Number=.,Type=Integer,Description=\"Spanning paired-read support for the ref and alt alleles in the order listed\">\n";
  os << "##FORMAT=<ID=SR,Number=.,Type=Integer,Description=\"Split reads for the ref and alt alleles in the order listed\">\n";
}

void VcfWriterRnaSV::addHeaderFilters(std::ostream& os) const
{
  os << "##FILTER=<ID=" << SVScoreInfoRna::rnaFilterLabel
     << ",Description=\"RNA fusion calls without both split read and spanning pair support\">\n";
  os << "##FILTER=<ID=" << SVScoreInfoRna::impreciseLabel
     << ",Description=\"RNA fusion candidates for which no spanning contig was found\">\n";
  os << "##FILTER=<ID=" << SVScoreInfoRna::localLabel
     << ",Description=\"RNA call covering short genomic distance\">\n";
}

void VcfWriterRnaSV::modifyTranslocInfo(
    const SVCandidate&             sv,
    const SVScoreInfo*             baseScoringInfoPtr,
    const bool                     isFirstOfPair,
    const SVCandidateAssemblyData& assemblyData,
    InfoTag_t&                     infotags) const
{
  assert(baseScoringInfoPtr);
  const SVScoreInfo& baseScoringInfo(*baseScoringInfoPtr);

  infotags.push_back(
      str(boost::format("BND_DEPTH=%i") %
          (isFirstOfPair ? baseScoringInfo.bp1MaxDepth : baseScoringInfo.bp2MaxDepth)));
  infotags.push_back(
      str(boost::format("MATE_BND_DEPTH=%i") %
          (isFirstOfPair ? baseScoringInfo.bp2MaxDepth : baseScoringInfo.bp1MaxDepth)));
  {
    /// TODO better multisample handler here:
    const unsigned sampleIndex(0);

    const SVSampleAlleleInfo& refinfo(baseScoringInfo.samples[sampleIndex].ref);
    infotags.push_back(
        str(boost::format("REF_COUNT=%i") % (isFirstOfPair ? refinfo.confidentSplitReadAndPairCountRefBp1
                                                           : refinfo.confidentSplitReadAndPairCountRefBp2)));
    infotags.push_back(str(
        boost::format("MATE_REF_COUNT=%i") % (isFirstOfPair ? refinfo.confidentSplitReadAndPairCountRefBp2
                                                            : refinfo.confidentSplitReadAndPairCountRefBp1)));
  }
  {
    // if (!assemblyData.isSpanning) return;

    const bool isFirst = (assemblyData.bporient.isBp1First == isFirstOfPair);
    if (isFirst) infotags.push_back("RNA_FIRST");
    if (assemblyData.bporient.isTranscriptStrandKnown) infotags.push_back("RNA_STRANDED");

    if (!isFirstOfPair)
      return;  // only the first breakpoint gets the additional RNA info attached to its VCF entry

    infotags.push_back(
        str(boost::format("RNA_FwRvReads=%i,%i") % sv.forwardTranscriptStrandReadCount %
            sv.reverseTranscriptStrandReadCount));
    infotags.push_back(str(boost::format("RNA_Reads=%i") % sv.bp2.lowresEvidence.getTotal()));
    const unsigned numContigs(assemblyData.contigs.size());
    if (numContigs > 0) {
      if (numContigs != assemblyData.spanningAlignments.size())
        infotags.push_back(
            str(boost::format("ERROR=%i,%i") % numContigs % assemblyData.spanningAlignments.size()));
      const unsigned int bestAlignmentIdx(assemblyData.bestAlignmentIndex);
      if (numContigs <= bestAlignmentIdx)
        infotags.push_back(str(boost::format("ERROR2=%i,%i") % numContigs % bestAlignmentIdx));
      infotags.push_back(str(boost::format("RNA_CONTIG=%s") % assemblyData.contigs[bestAlignmentIdx].seq));
      const auto& bestAlignment(assemblyData.spanningAlignments[bestAlignmentIdx]);
      infotags.push_back(
          str(boost::format("RNA_CONTIG_ALN=%i,%i") % apath_matched_length(bestAlignment.align1.apath) %
              apath_matched_length(bestAlignment.align2.apath)));
    }
  }
#ifdef DEBUG_VCF
  addDebugInfo(isFirstOfPair, sv, assemblyData, infotags);
#endif
}

#ifdef DEBUG_VCF
static void addDebugInfo(
    const bool                     isFirstOfPair,
    const SVCandidate&             sv,
    const SVCandidateAssemblyData& assemblyData,
    VcfWriterSV::InfoTag_t&        infotags)
{
  if (!assemblyData.isSpanning) return;

  const bool isFirst     = (assemblyData.bporient.isBp1First == isFirstOfPair);
  const bool isRightOpen = (isFirstOfPair ? sv.bp1.state : sv.bp2.state) == SVBreakendState::RIGHT_OPEN;
  infotags.push_back(str(boost::format("FOOBAR_FW=%1%") % (isFirst == isRightOpen)));

  if (!isFirst) return;  // only the first breakpoint gets the alignments attached to its VCF entry

  infotags.push_back(str(boost::format("FOOBAR_bp1=%i;bp2=%i") % sv.bp1.interval.tid % sv.bp2.interval.tid));

  // there can be several contigs per breakend, so we iterate over all of them.
  const unsigned numContigs(assemblyData.contigs.size());
  infotags.push_back(str(boost::format("FOOBAR_NCONTIGS=%i") % numContigs));
  if (numContigs > 0) {
    if (numContigs != assemblyData.spanningAlignments.size())
      infotags.push_back(
          str(boost::format("FOOBAR_ERROR=%i;%i") % numContigs % assemblyData.spanningAlignments.size()));
    if (numContigs <= assemblyData.bestAlignmentIndex)
      infotags.push_back(
          str(boost::format("FOOBAR_ERROR2=%i;%i") % numContigs % assemblyData.bestAlignmentIndex));

    infotags.push_back(str(boost::format("FOOBAR_BEST=%i") % assemblyData.bestAlignmentIndex));
    //infotags.push_back(str(boost::format("FOOBAR_EXTCONTIG=%s") % assemblyData.extendedContigs[assemblyData.bestAlignmentIndex]));
    infotags.push_back(
        str(boost::format("FOOBAR_CONTIGcount=%i") %
            assemblyData.contigs[assemblyData.bestAlignmentIndex].supportReads.size()));
  }
}
#endif

void VcfWriterRnaSV::writeFilter(const boost::any specializedScoringInfo, std::ostream& os) const
{
  const SVScoreInfoRna& rnaScoringInfo(*boost::any_cast<const SVScoreInfoRna*>(specializedScoringInfo));
  writeFilters(rnaScoringInfo.filters, os);
}

void VcfWriterRnaSV::modifySample(
    const SVCandidate& sv,
    const SVScoreInfo* baseScoringInfoPtr,
    const boost::any /*specializedScoringInfo*/,
    SampleTag_t& sampletags) const
{
  assert(baseScoringInfoPtr);
  const SVScoreInfo& baseScoringInfo(*baseScoringInfoPtr);
  const unsigned     sampleCount(baseScoringInfo.samples.size());

  std::vector<std::string> values(sampleCount);
  for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex) {
    const SVSampleInfo& sampleInfo(baseScoringInfo.samples[sampleIndex]);
    values[sampleIndex] =
        str(boost::format("%i,%i") % sampleInfo.ref.spanningPairCount % sampleInfo.alt.spanningPairCount);
  }
  sampletags.push_back(std::make_pair("PR", values));

  if (sv.isImprecise()) return;

  for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex) {
    const SVSampleInfo& sampleInfo(baseScoringInfo.samples[sampleIndex]);
    values[sampleIndex] =
        str(boost::format("%i,%i") % sampleInfo.ref.splitReadCount % sampleInfo.alt.splitReadCount);
  }
  sampletags.push_back(std::make_pair("SR", values));
}

void VcfWriterRnaSV::writeSV(
    const SVCandidateSetData&      svData,
    const SVCandidateAssemblyData& adata,
    const SVCandidate&             sv,
    const SVId&                    svId,
    const SVScoreInfo&             baseScoringInfo,
    const SVScoreInfoRna&          rnaInfo,
    const EventInfo&               event) const
{
  writeSVCore(svData, adata, sv, svId, &baseScoringInfo, &rnaInfo, event, true);
}
