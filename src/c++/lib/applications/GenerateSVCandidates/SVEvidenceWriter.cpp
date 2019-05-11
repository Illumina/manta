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
/// \author Xiaoyu Chen
///

#include "SVEvidenceWriter.hpp"

#include <iostream>

#include "manta/BamStreamerUtils.hpp"

//#define DEBUG_SUPPORT

#ifdef DEBUG_SUPPORT
#include "blt_util/log.hpp"
#endif

std::ostream& operator<<(std::ostream& os, const SVEvidenceWriterRead& read)
{
  os << read.tid << ":" << read.pos << "\t";
  for (const auto& sv : read.SVs) {
    os << sv.first << ",";
  }

  return os;
}

std::ostream& operator<<(std::ostream& os, const SVEvidenceWriterReadPair& readPair)
{
  os << readPair.read1 << "\n" << readPair.read2 << "\n";
  return os;
}

std::ostream& operator<<(std::ostream& os, const SVEvidenceWriterSampleData& sample)
{
  for (const auto& frg : sample.supportFrags) {
    os << "qname =" << frg.first << "\n" << frg.second;
  }

  return os;
}

std::ostream& operator<<(std::ostream& os, const SVEvidenceWriterData& data)
{
  unsigned size(data.sampleData.size());
  for (unsigned i = 0; i < size; i++) {
    os << "sample index = " << i << "\n" << data.sampleData[i];
  }

  return os;
}

void SVEvidenceWriter::processBamRecords(
    bam_streamer&              origBamStream,
    const GenomeInterval&      interval,
    const support_fragments_t& supportFrags,
    SynchronizedBamWriter&     bamWriter)
{
#ifdef DEBUG_SUPPORT
  log_os << __FUNCTION__ << "  target interval: " << interval << "\n";
#endif

  origBamStream.resetRegion(interval.tid, interval.range.begin_pos(), interval.range.end_pos());
  while (origBamStream.next()) {
    const bam_record* origBamRec(origBamStream.get_record_ptr());
    bam_record        bamRec(*origBamRec);

    const std::string                   qname(bamRec.qname());
    support_fragments_t::const_iterator suppFragsIter(supportFrags.find(qname));
    if (suppFragsIter != supportFrags.end()) {
      const SVEvidenceWriterReadPair& supportFrag(suppFragsIter->second);
      const bool                      isR1Matched(
          bamRec.is_first() && (bamRec.target_id() == supportFrag.read1.tid) &&
          (bamRec.pos() == supportFrag.read1.pos));
      const bool isR2Matched(
          (!bamRec.is_first()) && (bamRec.target_id() == supportFrag.read2.tid) &&
          (bamRec.pos() == supportFrag.read2.pos));

      if (isR1Matched || isR2Matched) {
        const SVEvidenceWriterRead& read(isR1Matched ? supportFrag.read1 : supportFrag.read2);
#ifdef DEBUG_SUPPORT
        log_os << __FUNCTION__ << "  matched supporting read: " << read << "\n";
#endif

        bam1_t& br(*(bamRec.get_data()));
        // add new customized field of SV IDs that the read supports
        bool        isFirst(true);
        std::string svStr;
        for (const auto& sv : read.SVs) {
          if (!isFirst) svStr.append(",");
          svStr.append(sv.first);
          for (const auto& svType : sv.second) {
            svStr.append('|' + svType);
          }
          if (isFirst) isFirst = false;
        }

        static const char svtag[] = {'Z', 'M'};
        // Do lots of ugly casting on svStr to fit htsapi signature. Usage is actually const in htslib:
        bam_aux_append(
            &br,
            svtag,
            'Z',
            (svStr.size() + 1),
            reinterpret_cast<uint8_t*>(const_cast<char*>(svStr.c_str())));

        // Update bam record bin value
        bam_update_bin(br);
        // write to bam

        bamWriter.put_record(&br);
      }
    }
  }
}

void SVEvidenceWriter::writeSupportBam(
    const bam_streamer_ptr&           origBamStreamPtr,
    const SVEvidenceWriterSampleData& svSupportFrags,
    SynchronizedBamWriter&            bamWriter)
{
  std::vector<SVEvidenceWriterRead> supportReads;
  const support_fragments_t&        supportFrags(svSupportFrags.supportFrags);
  for (const auto& frg : supportFrags) {
    supportReads.push_back(frg.second.read1);
    supportReads.push_back(frg.second.read2);
  }
  // sort all the reads w.r.t. genomic positions
  std::sort(supportReads.begin(), supportReads.end());

  // generate a set of intervals containing overlapping reads
  const int                   readDistance(100);
  int                         lastTid = -1;
  int                         lastPos = -1;
  std::vector<GenomeInterval> intervals;
  for (const auto& suppRd : supportReads) {
    if ((lastTid == suppRd.tid) && (lastPos + readDistance >= suppRd.pos)) {
      GenomeInterval& interval(intervals.back());
      interval.range.set_end_pos(suppRd.pos);
    } else {
      GenomeInterval interval(suppRd.tid, suppRd.pos - 1, suppRd.pos);
      intervals.push_back(interval);
    }

    lastTid = suppRd.tid;
    lastPos = suppRd.pos;
  }

  bam_streamer& origBamStream(*origBamStreamPtr);
  for (const auto& interval : intervals) {
    processBamRecords(origBamStream, interval, supportFrags, bamWriter);
  }
}

SVEvidenceWriterSharedData::SVEvidenceWriterSharedData(const GSCOptions& opt)
{
  if (!opt.isGenerateEvidenceBam()) return;

  const unsigned sampleSize(opt.alignFileOpt.alignmentFilenames.size());
  for (unsigned sampleIndex(0); sampleIndex < sampleSize; ++sampleIndex) {
    const std::string  evidenceBamName(opt.evidenceBamStub + ".bam_" + std::to_string(sampleIndex) + ".bam");
    const bam_streamer bamStreamer(
        opt.alignFileOpt.alignmentFilenames[sampleIndex].c_str(), opt.referenceFilename.c_str());

    m_evidenceBamWriterPtrs.push_back(
        std::make_shared<SynchronizedBamWriter>(evidenceBamName.c_str(), bamStreamer.get_header()));
  }
}

SVEvidenceWriter::SVEvidenceWriter(
    const GSCOptions& opt, std::shared_ptr<SVEvidenceWriterSharedData> sharedData)
  : m_isGenerateEvidenceBam(opt.isGenerateEvidenceBam()),
    m_sampleSize(opt.alignFileOpt.alignmentFilenames.size()),
    m_sharedData(sharedData)
{
  if (!m_isGenerateEvidenceBam) return;

  openBamStreams(opt.referenceFilename, opt.alignFileOpt.alignmentFilenames, m_origBamStreamPtrs);
  assert(m_origBamStreamPtrs.size() == m_sampleSize);
}

void SVEvidenceWriter::write(const SVEvidenceWriterData& svEvidenceWriterData)
{
  if (!m_isGenerateEvidenceBam) return;

  for (unsigned sampleIndex(0); sampleIndex < m_sampleSize; ++sampleIndex) {
    writeSupportBam(
        m_origBamStreamPtrs[sampleIndex],
        svEvidenceWriterData.sampleData[sampleIndex],
        m_sharedData->getBamWriter(sampleIndex));
  }
}
