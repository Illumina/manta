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
/// \author Chris Saunders
/// \author Felix Schlesinger
/// \author Naoki Nariai
/// \author Xiaoyu Chen
///

#include "format/VcfWriterSV.hpp"

#include "blt_util/log.hpp"
#include "blt_util/seq_util.hpp"
#include "blt_util/string_util.hpp"
#include "common/Exceptions.hpp"
#include "htsapi/samtools_fasta_util.hpp"
#include "htsapi/vcf_util.hpp"
#include "manta/SVCandidateUtil.hpp"

#include <iostream>
#include <sstream>

//#define DEBUG_VCF

#ifdef DEBUG_VCF
#include "blt_util/log.hpp"
#endif

VcfWriterSV::VcfWriterSV(
    const std::string&     referenceFilename,
    const bam_header_info& bamHeaderInfo,
    const std::string&     outputFilename,
    const bool&            isOutputContig)
  : _referenceFilename(referenceFilename),
    _isOutputContig(isOutputContig),
    _stream(outputFilename),
    _header(bamHeaderInfo)
{
}

void VcfWriterSV::writeHeader(
    const char* progName, const char* progVersion, const std::vector<std::string>& sampleNames) const
{
  std::ostringstream oss;
  writeHeaderPrefix(progName, progVersion, oss);
  writeHeaderColumnKey(sampleNames, oss);
  _stream.write(oss.str());
}

void VcfWriterSV::writeHeaderPrefix(const char* progName, const char* progVersion, std::ostream& os) const
{
  os << "##fileformat=VCFv4.1\n";
  os << "##fileDate=" << vcf_fileDate << "\n";
  os << "##source=" << progName << " " << progVersion << "\n";
  os << "##reference=file://" << _referenceFilename << "\n";

  for (const bam_header_info::chrom_info& cdata : _header.chrom_data) {
    os << "##contig=<ID=" << cdata.label << ",length=" << cdata.length << ">\n";
  }

  /// vcf 4.1 reserved/suggested INFO tags:
  os << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n";
  os << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
  os << "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n";
  os << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n";
  os << "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS\">\n";
  os << "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END\">\n";
  os << "##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"CIGAR alignment for each alternate indel allele\">\n";
  os << "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakend\">\n";
  os << "##INFO=<ID=EVENT,Number=1,Type=String,Description=\"ID of event associated to breakend\">\n";
  os << "##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical homology at event breakpoints\">\n";
  os << "##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical homology at event breakpoints\">\n";

  /// custom INFO tags:
  os << "##INFO=<ID=SVINSLEN,Number=.,Type=Integer,Description=\"Length of insertion\">\n";
  os << "##INFO=<ID=SVINSSEQ,Number=.,Type=String,Description=\"Sequence of insertion\">\n";
  os << "##INFO=<ID=LEFT_SVINSSEQ,Number=.,Type=String,Description=\"Known left side of insertion for an insertion of unknown length\">\n";
  os << "##INFO=<ID=RIGHT_SVINSSEQ,Number=.,Type=String,Description=\"Known right side of insertion for an insertion of unknown length\">\n";

  // if "--outputContig" is specified, then print out INFO tag for Assembled contig sequence
  if (_isOutputContig) {
    os << "##INFO=<ID=CONTIG,Number=1,Type=String,Description=\"Assembled contig sequence\">\n";
  }

  addHeaderInfo(os);

  addHeaderFormat(os);

  addHeaderFilters(os);

  os << "##ALT=<ID=DEL,Description=\"Deletion\">\n";
  os << "##ALT=<ID=INS,Description=\"Insertion\">\n";
  os << "##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">\n";
}

static void writeHeaderColKeyPrefix(std::ostream& os)
{
  os << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
}

void VcfWriterSV::writeHeaderColumnKey(const std::vector<std::string>& sampleNames, std::ostream& os) const
{
  writeHeaderColKeyPrefix(os);
  if (!sampleNames.empty()) {
    os << "\tFORMAT";

    for (const std::string& sampleName : sampleNames) {
      os << '\t' << sampleName;
    }
  }
  os << '\n';
}

static void makeInfoField(const VcfWriterSV::InfoTag_t& info, std::ostream& os)
{
  static const char sep(';');
  bool              isFirst(true);
  for (const std::string& is : info) {
    if (!isFirst)
      os << sep;
    else
      isFirst = false;
    os << is;
  }
}

static void makeFormatSampleField(const VcfWriterSV::SampleTag_t& sample, std::ostream& os)
{
  static const char sep(':');

  if (sample.empty()) return;

  {
    // first write FORMAT field:
    os << '\t';

    bool isFirst(true);
    for (const VcfWriterSV::SampleTag_t::value_type& fs : sample) {
      if (!isFirst)
        os << sep;
      else
        isFirst = false;

      assert(!fs.first.empty());
      os << fs.first;
    }
  }

  unsigned nSamples(0);
  for (const VcfWriterSV::SampleTag_t::value_type& fs : sample) {
    const unsigned ns(fs.second.size());
    nSamples = std::max(nSamples, ns);
  }

  for (unsigned sampleIndex(0); sampleIndex < nSamples; ++sampleIndex) {
    os << '\t';

    // next write SAMPLE field:
    {
      bool isFirst(true);
      for (const VcfWriterSV::SampleTag_t::value_type& fs : sample) {
        if (!isFirst)
          os << sep;
        else
          isFirst = false;

        if (fs.second.size() <= sampleIndex) {
          os << '.';
        } else if (fs.second[sampleIndex].empty()) {
          os << '.';
        } else {
          os << fs.second[sampleIndex];
        }
      }
    }
  }
}

#ifdef DEBUG_VCF
static void addDebugInfo(
    const SVBreakend&              bp1,
    const SVBreakend&              bp2,
    const bool                     isFirstOfPair,
    const SVCandidateAssemblyData& assemblyData,
    VcfWriterSV::InfoTag_t&        infotags)
{
  if (!isFirstOfPair) return;

  // store alignment start + cigar string for each section of the jumping alignment.
  // there can be several contigs per breakend, so we iterate over all of them.
  // only the first breakpoint gets the alignments attached to its VCF entry

  if (assemblyData.isSpanning) {
    const unsigned numAlign(assemblyData.spanningAlignments.size());
    std::string    cigar1;
    std::string    cigar2;
    for (unsigned alignIndex(0); alignIndex < numAlign; ++alignIndex) {
      const SVCandidateAssemblyData::JumpAlignmentResultType align(
          assemblyData.spanningAlignments[alignIndex]);
      infotags.push_back(
          str(boost::format("CTG_JALIGN_%i_POS_A=%d") % alignIndex %
              (bp1.interval.range.begin_pos() + align.align1.beginPos)));
      infotags.push_back(
          str(boost::format("CTG_JALIGN_%i_POS_B=%d") % alignIndex %
              (bp2.interval.range.begin_pos() + align.align2.beginPos)));

      apath_to_cigar(align.align1.apath, cigar1);
      apath_to_cigar(align.align2.apath, cigar2);

      infotags.push_back(str(boost::format("CTG_JALIGN_%i_CIGAR_A=%s") % alignIndex % cigar1));
      infotags.push_back(str(boost::format("CTG_JALIGN_%i_CIGAR_B=%s") % alignIndex % cigar2));
    }
    const unsigned numContigs(assemblyData.contigs.size());

    infotags.push_back(str(boost::format("DEBUG_NCONTIGS=%i") % numContigs));
    infotags.push_back(str(
        boost::format("DEBUG_BESTContig=%s") % assemblyData.contigs[assemblyData.bestAlignmentIndex].seq));
    infotags.push_back(
        str(boost::format("DEBUG_CONTIGReads=%i") %
            assemblyData.contigs[assemblyData.bestAlignmentIndex].supportReads.size()));
    infotags.push_back(str(
        boost::format("DEBUG_CONTIGLeftAln=%i") %
        apath_matched_length(assemblyData.spanningAlignments[assemblyData.bestAlignmentIndex].align1.apath)));
    infotags.push_back(str(
        boost::format("DEBUG_CONTIGRightAln=%i") %
        apath_matched_length(assemblyData.spanningAlignments[assemblyData.bestAlignmentIndex].align2.apath)));
  }
}
#endif

static void addSharedInfo(const EventInfo& event, VcfWriterSV::InfoTag_t& infoTags)
{
  if (event.isEvent()) {
    infoTags.push_back(str(boost::format("EVENT=%i") % event.label));
  }
}

/// Add HOMLEN & HOMSEQ info tags for breakend homology
///
/// \param[in] bpRange breakpoint homology region
/// \param[in] bpPosAdjust breakpoint position adjustment according to breakend direction to match vcf spec
/// \param infoTags append HOMELEN & HOMSEQ to infoTags
static void addHomologyInfo(
    const std::string&      refFile,
    const std::string&      chrom,
    const known_pos_range2& bpRange,
    const pos_t             bpPosAdjust,
    VcfWriterSV::InfoTag_t& infoTags)
{
  if (bpRange.size() > 1) {
    // get breakend homology sequence
    infoTags.push_back(str(boost::format("HOMLEN=%i") % (bpRange.size() - 1)));
    std::string homSeq;
    // the get region seq function below takes closed-closed, 0-based endpoints
    const int homBegin(bpRange.begin_pos() + bpPosAdjust + 1);
    const int homEnd(bpRange.end_pos() + bpPosAdjust - 1);
    get_standardized_region_seq(refFile, chrom, homBegin, homEnd, homSeq);
    infoTags.push_back(str(boost::format("HOMSEQ=%s") % (homSeq)));
  }
}

void VcfWriterSV::writeTransloc(
    const SVCandidate& sv,
    const SVId&        svId,
    const SVScoreInfo* baseScoringInfoPtr,
    const boost::any   specializedScoringInfo,
    const bool         isFirstBreakend,
    const SVCandidateSetData& /*svData*/,
    const SVCandidateAssemblyData& adata,
    const EventInfo&               event) const
{
  const bool isImprecise(sv.isImprecise());
  const bool isBreakendRangeSameShift(sv.isBreakendRangeSameShift());

  const SVBreakend& bpA(isFirstBreakend ? sv.bp1 : sv.bp2);
  const SVBreakend& bpB(isFirstBreakend ? sv.bp2 : sv.bp1);

  InfoTag_t   infoTags;
  SampleTag_t sampleTags;

  // get CHROM
  const std::string& chrom(_header.chrom_data[bpA.interval.tid].label);
  const std::string& mateChrom(_header.chrom_data[bpB.interval.tid].label);

  const known_pos_range2& bpARange(bpA.interval.range);
  const known_pos_range2& bpBRange(bpB.interval.range);

  if (!isImprecise) {
    assert(bpARange.size() == bpBRange.size());
  }

  // get POS
  pos_t pos(bpARange.center_pos() + 1);
  pos_t matePos(bpBRange.center_pos() + 1);
  if (!isImprecise) {
    pos = bpARange.begin_pos() + 1;
    if (isBreakendRangeSameShift) {
      matePos = bpBRange.begin_pos() + 1;
    } else {
      matePos = bpBRange.end_pos();
    }
  }

  // TODO: improve circular genome handler:
  if ((pos < 1) || (matePos < 1)) return;

  // get ID
  const std::string& localId(isFirstBreakend ? svId.localId : svId.mateId);
  const std::string& mateId(isFirstBreakend ? svId.mateId : svId.localId);

  // get REF
  std::string ref;
  get_standardized_region_seq(_referenceFilename, chrom, pos - 1, pos - 1, ref);

  assert(1 == ref.size());

  const bool         isReverseInsertSeq(!(isFirstBreakend || (bpA.state != bpB.state)));
  std::string        tmpString;
  const std::string* insertSeqPtr(&sv.insertSeq);
  if (isReverseInsertSeq) {
    tmpString    = reverseCompCopyStr(sv.insertSeq);
    insertSeqPtr = &tmpString;
  }
  const std::string& insertSeq(*insertSeqPtr);

  // build alt:
  boost::format altFormat("%4%%3%%1%:%2%%3%%5%");
  {
    std::string altPrefix;
    std::string altSuffix;
    if (bpA.state == SVBreakendState::RIGHT_OPEN) {
      altPrefix = ref + insertSeq;
    } else if (bpA.state == SVBreakendState::LEFT_OPEN) {
      altSuffix = insertSeq + ref;
    } else {
      assert(false && "Unexpected bpA.state");
    }

    char altSep('?');
    if (bpB.state == SVBreakendState::RIGHT_OPEN) {
      altSep = ']';
    } else if (bpB.state == SVBreakendState::LEFT_OPEN) {
      altSep = '[';
    } else {
      assert(false && "Unexpected bpB.state");
    }

    altFormat % mateChrom % matePos % altSep % altPrefix % altSuffix;
  }

  // build INFO field
  infoTags.push_back("SVTYPE=BND");
  infoTags.push_back("MATEID=" + mateId);
  if (isImprecise) {
    infoTags.push_back("IMPRECISE");
  } else if (_isOutputContig) {
    infoTags.push_back("CONTIG=" + sv.contigSeq);
  }

  if (bpARange.size() > 1) {
    infoTags.push_back(
        str(boost::format("CIPOS=%i,%i") % ((bpARange.begin_pos() + 1) - pos) % (bpARange.end_pos() - pos)));
  }

  if (!isImprecise) {
    const pos_t bpAPosAdjust(0);
    addHomologyInfo(_referenceFilename, chrom, bpARange, bpAPosAdjust, infoTags);
  }

  if (!insertSeq.empty()) {
    infoTags.push_back(str(boost::format("SVINSLEN=%i") % (insertSeq.size())));
    infoTags.push_back(str(boost::format("SVINSSEQ=%s") % (insertSeq)));
  }

  addSharedInfo(event, infoTags);

  modifyInfo(event, specializedScoringInfo, infoTags);
  modifyTranslocInfo(sv, baseScoringInfoPtr, isFirstBreakend, adata, infoTags);

  modifySample(sv, baseScoringInfoPtr, specializedScoringInfo, sampleTags);
#ifdef DEBUG_VCF
  addDebugInfo(bpA, bpB, isFirstBreakend, adata, infoTags);
#endif

  // write out record:
  std::ostringstream oss;
  oss << chrom << '\t' << pos << '\t' << localId  // ID
      << '\t' << ref                              // REF
      << '\t' << str(altFormat)                   // ALT
      << '\t';
  writeQual(specializedScoringInfo, oss);
  oss << '\t';
  writeFilter(specializedScoringInfo, oss);
  oss << '\t';
  makeInfoField(infoTags, oss);            // INFO
  makeFormatSampleField(sampleTags, oss);  // FORMAT + SAMPLE
  oss << '\n';
  _stream.write(oss.str());
}

void VcfWriterSV::writeTranslocPair(
    const SVCandidate&             sv,
    const SVId&                    svId,
    const SVScoreInfo*             baseScoringInfoPtr,
    const boost::any               specializedScoringInfo,
    const SVCandidateSetData&      svData,
    const SVCandidateAssemblyData& adata,
    const EventInfo&               event) const
{
  writeTransloc(sv, svId, baseScoringInfoPtr, specializedScoringInfo, true, svData, adata, event);
  writeTransloc(sv, svId, baseScoringInfoPtr, specializedScoringInfo, false, svData, adata, event);
}

void VcfWriterSV::writeIndel(
    const SVCandidate& sv,
    const SVId&        svId,
    const SVScoreInfo* baseScoringInfoPtr,
    const boost::any   specializedScoringInfo,
    const bool         isIndel,
    const EventInfo&   event) const
{
  const bool isImprecise(sv.isImprecise());
  const bool isBreakendRangeSameShift(sv.isBreakendRangeSameShift());

  const bool isBp1First(sv.bp1.interval.range.begin_pos() <= sv.bp2.interval.range.begin_pos());

  const SVBreakend& bpA(isBp1First ? sv.bp1 : sv.bp2);
  const SVBreakend& bpB(isBp1First ? sv.bp2 : sv.bp1);

  InfoTag_t   infoTags;
  SampleTag_t sampleTags;

  // get CHROM
  const std::string& chrom(_header.chrom_data[sv.bp1.interval.tid].label);

  const known_pos_range2& bpARange(bpA.interval.range);
  const known_pos_range2& bpBRange(bpB.interval.range);

  if (!isImprecise) {
    assert(bpARange.size() == bpBRange.size());
  }

  // above this size all records use symbolic alleles (ie. <DEL>):
  static const unsigned maxNonSymbolicRecordSize(1000);

  // if the variant is a combination of simple insertion and deletions, and below
  // a large-event size threshold, it is classified as a small variant. In this case
  // we report the event using full REF and ALT sequences, plus a CIGAR string for
  // complex in/del combinations
  //
  bool isSmallVariant(false);
  if ((!isImprecise) && isIndel && (!sv.isUnknownSizeInsertion)) {
    const unsigned deleteSize(bpBRange.begin_pos() - bpARange.begin_pos());
    const unsigned insertSize(sv.insertSeq.size());

    const bool isSmallDelete(deleteSize <= maxNonSymbolicRecordSize);
    const bool isSmallInsert(insertSize <= maxNonSymbolicRecordSize);

    isSmallVariant = (isSmallDelete && isSmallInsert);
  }

  // get POS and endPos,
  // first compute internal coordinates and then transform per vcf conventions:
  pos_t internal_pos(bpARange.center_pos());
  pos_t internal_endPos(bpBRange.center_pos());
  if (!isImprecise) {
    internal_pos = bpARange.begin_pos();
    if (isBreakendRangeSameShift) {
      internal_endPos = bpBRange.begin_pos();
    } else {
      internal_endPos = (bpBRange.end_pos() - 1);
    }
  }

  // now create external pos values for vcf only
  // everything is +1'd to get out zero-indexed coordinates:
  pos_t pos(internal_pos + 1);
  pos_t endPos(internal_endPos + 1);

  // variants are adjusted by up to one base according to breakend direction to match vcf spec:
  const pos_t bpAPosAdjust(bpA.getLeftSideOfBkptAdjustment());
  const pos_t bpBPosAdjust(bpB.getLeftSideOfBkptAdjustment());
  pos += bpAPosAdjust;
  endPos += bpBPosAdjust;

  if (pos < 1) return;

  // get REF
  std::string ref;
  {
    const pos_t beginRefPos(pos - 1);
    pos_t       endRefPos(beginRefPos);
    if (isSmallVariant) endRefPos = endPos - 1;

    get_standardized_region_seq(_referenceFilename, chrom, beginRefPos, endRefPos, ref);

    if (static_cast<unsigned>(1 + endRefPos - beginRefPos) != ref.size()) {
      using namespace illumina::common;

      std::ostringstream oss;
      oss << "Unexpected reference allele size: " << ref.size() << "\n";
      oss << "\tExpected: " << (1 + endRefPos - beginRefPos) << "\n";
      oss << "\tbeginRefPos: " << beginRefPos << " endRefPos: " << endRefPos
          << " isSmallVariant: " << isSmallVariant << "\n";
      BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
    }
  }

  // build alt:
  std::string alt;
  if (isSmallVariant) {
    alt = ref[0] + sv.insertSeq;
  } else {
    alt = str(boost::format("<%s>") % svId.getLabel());
  }

  // build INFO field
  std::vector<std::string> words;
  split_string(svId.getLabel(), ':', words);
  {
    // note that there's a reasonable argument for displaying these tags only when a
    // symbolic allele is used (by a strict reading of the vcf spec) -- we instead
    // print these fields for all variants for uniformity within the manta vcf:
    //
    infoTags.push_back(str(boost::format("END=%i") % endPos));
    infoTags.push_back(str(boost::format("SVTYPE=%s") % words[0]));
    const pos_t refLen(endPos - pos);
    pos_t       svLen(refLen);

    if (!sv.isUnknownSizeInsertion) {
      if (isIndel) {
        const pos_t insertLen(static_cast<pos_t>(sv.insertSeq.size()));
        if (insertLen > refLen) {
          svLen = insertLen;
        } else {
          svLen = -refLen;
        }
      }
      infoTags.push_back(str(boost::format("SVLEN=%i") % (svLen)));
    }
  }

  if (isSmallVariant) {
    if (!sv.insertAlignment.empty()) {
      std::string cigar;
      apath_to_cigar(sv.insertAlignment, cigar);

      // add the 1M to signify the leading reference base:
      infoTags.push_back(str(boost::format("CIGAR=1M%s") % cigar));
    }
  }

  if (isImprecise) {
    infoTags.push_back("IMPRECISE");
  } else if (_isOutputContig) {
    infoTags.push_back("CONTIG=" + sv.contigSeq);
  }

  if (bpARange.size() > 1) {
    infoTags.push_back(
        str(boost::format("CIPOS=%i,%i") % (bpARange.begin_pos() - internal_pos) %
            ((bpARange.end_pos() - 1) - internal_pos)));
  }

  if (!isSmallVariant) {
    if (bpBRange.size() > 1) {
      infoTags.push_back(
          str(boost::format("CIEND=%i,%i") % (bpBRange.begin_pos() - internal_endPos) %
              ((bpBRange.end_pos() - 1) - internal_endPos)));
    }
  }

  if (!isImprecise) {
    addHomologyInfo(_referenceFilename, chrom, bpARange, bpAPosAdjust, infoTags);
  }

  if (!isSmallVariant) {
    if (!(sv.insertSeq.empty() || sv.isUnknownSizeInsertion)) {
      infoTags.push_back(str(boost::format("SVINSLEN=%i") % (sv.insertSeq.size())));
      if (isBp1First || (bpA.state != bpB.state)) {
        infoTags.push_back(str(boost::format("SVINSSEQ=%s") % (sv.insertSeq)));
      } else {
        infoTags.push_back(str(boost::format("SVINSSEQ=%s") % reverseCompCopyStr(sv.insertSeq)));
      }
    }
  }

  if (sv.isUnknownSizeInsertion) {
    if (!sv.unknownSizeInsertionLeftSeq.empty()) {
      infoTags.push_back(str(boost::format("LEFT_SVINSSEQ=%s") % (sv.unknownSizeInsertionLeftSeq)));
    }

    if (!sv.unknownSizeInsertionRightSeq.empty()) {
      infoTags.push_back(str(boost::format("RIGHT_SVINSSEQ=%s") % (sv.unknownSizeInsertionRightSeq)));
    }
  }

  addSharedInfo(event, infoTags);

  modifyInfo(event, specializedScoringInfo, infoTags);
  modifyInvdelInfo(sv, isBp1First, infoTags);

  modifySample(sv, baseScoringInfoPtr, specializedScoringInfo, sampleTags);

  // write out record:
  std::ostringstream oss;
  oss << chrom << '\t' << pos << '\t' << svId.localId  // ID
      << '\t' << ref                                   // REF
      << '\t' << alt                                   // ALT
      << '\t';
  writeQual(specializedScoringInfo, oss);
  oss << '\t';
  writeFilter(specializedScoringInfo, oss);
  oss << '\t';
  makeInfoField(infoTags, oss);            // INFO
  makeFormatSampleField(sampleTags, oss);  // FORMAT + SAMPLE
  oss << '\n';
  _stream.write(oss.str());
}

static bool isAcceptedSVType(const EXTENDED_SV_TYPE::index_t svType)
{
  using namespace EXTENDED_SV_TYPE;

  switch (svType) {
  case INTERTRANSLOC:
  case INTRATRANSLOC:
  case INVERSION:
  case INSERT:
  case DELETE:
  case TANDUP:
    return true;
  default:
    return false;
  }
}

void VcfWriterSV::writeSVCore(
    const SVCandidateSetData&      svData,
    const SVCandidateAssemblyData& adata,
    const SVCandidate&             sv,
    const SVId&                    svId,
    const SVScoreInfo*             baseScoringInfoPtr,
    const boost::any               specializedScoringInfo,
    const EventInfo&               event,
    const bool                     isForceIntraChromBnd) const
{
  using namespace EXTENDED_SV_TYPE;
  const index_t svType(getExtendedSVType(sv, isForceIntraChromBnd));

#ifdef DEBUG_VCF
  log_os << "VcfWriterSV::writeSVCore svType: " << EXTENDED_SV_TYPE::label(svType) << "\n";
#endif

  if (!isAcceptedSVType(svType)) {
    using namespace illumina::common;

    std::ostringstream oss;
    oss << "SV candidate cannot be classified: " << sv;
    BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
  }

  try {
    if (isSVTransloc(svType) || isSVInv(svType)) {
      writeTranslocPair(sv, svId, baseScoringInfoPtr, specializedScoringInfo, svData, adata, event);
    } else {
      const bool isIndel(isSVIndel(svType));
      writeIndel(sv, svId, baseScoringInfoPtr, specializedScoringInfo, isIndel, event);
    }
  } catch (...) {
    log_os << "Exception caught while attempting to write sv candidate to vcf: " << sv << "\n";
    log_os << "\tsvId: " << svId.getLabel() << " ext-svType: " << EXTENDED_SV_TYPE::label(svType) << "\n";
    throw;
  }
}

void VcfWriterSV::writeFilters(const std::set<std::string>& filters, std::ostream& os)
{
  if (filters.empty()) {
    os << "PASS";
  } else {
    bool isFirst(true);
    for (const std::string& filter : filters) {
      if (isFirst) {
        isFirst = false;
      } else {
        os << ';';
      }
      os << filter;
    }
  }
}

void VcfWriterSV::writeFilters(const std::set<std::string>& filters, std::string& s)
{
  std::ostringstream oss;
  writeFilters(filters, oss);
  s = oss.str();
}
