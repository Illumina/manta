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

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "GSCOptions.hpp"
#include "SynchronizedBamWriter.hpp"
#include "htsapi/bam_streamer.hpp"
#include "manta/SVCandidateSetData.hpp"

/// Records a single read that supports one or more SVs for evidence-BAM output
struct SVEvidenceWriterRead {
  typedef std::map<std::string, std::set<std::string>> SV_supportType_t;

  void addNewSV(const std::string& svId, const std::string& supportType)
  {
    if (SVs.find(svId) == SVs.end()) {
      std::set<std::string> supportTypeSet;
      SVs[svId] = supportTypeSet;
    }

    SVs[svId].insert(supportType);
  }

  bool operator<(const SVEvidenceWriterRead& rhs) const
  {
    if (tid < rhs.tid) return true;
    if (tid == rhs.tid) {
      return (pos < rhs.pos);
    }
    return false;
  }

  int              tid = -1;
  int              pos = 0;
  SV_supportType_t SVs;
};

std::ostream& operator<<(std::ostream& os, const SVEvidenceWriterRead& read);

/// Records a single fragment (read1 & read2)
/// that supports one or more SVs for evidence-BAM output,
/// indicating the evidence type
struct SVEvidenceWriterReadPair {
  void setReads(const bam_record& bamRead)
  {
    if (bamRead.is_first()) {
      read1.tid = bamRead.target_id();
      read1.pos = bamRead.pos();
      read2.tid = bamRead.mate_target_id();
      read2.pos = bamRead.mate_pos();
    } else if (bamRead.is_second()) {
      read1.tid = bamRead.mate_target_id();
      read1.pos = bamRead.mate_pos();
      read2.tid = bamRead.target_id();
      read2.pos = bamRead.pos();
    }
  }

  void addSpanningSupport(const std::string& svID)
  {
    read1.addNewSV(svID, "PR");
    read2.addNewSV(svID, "PR");
  }

  void addSplitSupport(const bool isRead1, const std::string& svID)
  {
    if (isRead1) {
      read1.addNewSV(svID, "SR");
      read2.addNewSV(svID, "SRM");
    } else {
      read2.addNewSV(svID, "SR");
      read1.addNewSV(svID, "SRM");
    }
  }

  SVEvidenceWriterRead read1;
  SVEvidenceWriterRead read2;
};

std::ostream& operator<<(std::ostream& os, const SVEvidenceWriterReadPair& readPair);

typedef std::map<std::string, SVEvidenceWriterReadPair> support_fragments_t;

/// Records all supporting fragments
/// that support one or more SVs for evidence-BAM output
struct SVEvidenceWriterSampleData {
  void clear() { supportFrags.clear(); }

  SVEvidenceWriterReadPair& getSupportFragment(const bam_record& bamRead)
  {
    const std::string qname(bamRead.qname());

    // create a new entry in the map
    if (supportFrags.find(qname) == supportFrags.end()) {
      SVEvidenceWriterReadPair newFrag;
      newFrag.setReads(bamRead);
      supportFrags[qname] = newFrag;
    }

    return supportFrags[qname];
  }

  SVEvidenceWriterReadPair& getSupportFragment(const SVCandidateSetSequenceFragment& seqFrag)
  {
    // Tentatively add an assertion
    // \TODO add the logic if only supplementary (NOT primary) reads being set
    assert(seqFrag.read1.isSet() || seqFrag.read2.isSet());

    const SVCandidateSetRead& read(seqFrag.read1.isSet() ? seqFrag.read1 : seqFrag.read2);

    const std::string qname(read.bamrec.qname());

    // create a new entry in the map
    if (supportFrags.find(qname) == supportFrags.end()) {
      SVEvidenceWriterReadPair newFrag;
      newFrag.setReads(read.bamrec);
      supportFrags[qname] = newFrag;
    }

    return supportFrags[qname];
  }

  support_fragments_t supportFrags;
};

std::ostream& operator<<(std::ostream& os, const SVEvidenceWriterSampleData& sample);

/// Records SV evidence fragments for all samples
struct SVEvidenceWriterData {
  explicit SVEvidenceWriterData(const unsigned sampleSize) : sampleData(sampleSize) {}

  void clear()
  {
    for (auto& data : sampleData) {
      data.clear();
    }
  }

  SVEvidenceWriterSampleData& getSampleData(const unsigned index)
  {
    assert(index < sampleData.size());
    return sampleData[index];
  }

  std::vector<SVEvidenceWriterSampleData> sampleData;
};

std::ostream& operator<<(std::ostream& os, const SVEvidenceWriterData& data);

/// Data that are shared by multiple SVEvidenceWriters (potentially operating on multiple threads)
class SVEvidenceWriterSharedData {
public:
  explicit SVEvidenceWriterSharedData(const GSCOptions& opt);

  SynchronizedBamWriter& getBamWriter(const unsigned bamIndex)
  {
    assert(bamIndex < m_evidenceBamWriterPtrs.size());
    assert(m_evidenceBamWriterPtrs[bamIndex]);
    return *(m_evidenceBamWriterPtrs[bamIndex]);
  }

private:
  std::vector<std::shared_ptr<SynchronizedBamWriter>> m_evidenceBamWriterPtrs;
};

/// \brief Coordinate all bookkeeping and data structures required to output evidence BAMs
///
/// Evidence BAMs are intended for visualization in debugging scenarios only. They are not part of the
/// standard calling workflow.
class SVEvidenceWriter {
public:
  SVEvidenceWriter(const GSCOptions& opt, std::shared_ptr<SVEvidenceWriterSharedData> sharedData);

  /// Write SV evidence reads into bam files
  void write(const SVEvidenceWriterData& svEvidenceWriterData);

  // Everything below is conceptually private, but kept here for unit testing until a fix can be written:

  typedef std::shared_ptr<bam_streamer> bam_streamer_ptr;

  static void processBamRecords(
      bam_streamer&              origBamStream,
      const GenomeInterval&      interval,
      const support_fragments_t& supportFrags,
      SynchronizedBamWriter&     bamWriter);

  static void writeSupportBam(
      const bam_streamer_ptr&           origBamStream,
      const SVEvidenceWriterSampleData& svSupportFrags,
      SynchronizedBamWriter&            bamWriter);

private:
  bool                                        m_isGenerateEvidenceBam;
  unsigned                                    m_sampleSize;
  std::vector<bam_streamer_ptr>               m_origBamStreamPtrs;
  std::shared_ptr<SVEvidenceWriterSharedData> m_sharedData;
};
