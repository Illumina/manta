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

#pragma once

#include "htsapi/bam_record.hpp"

#include <string>

/// this is the first of 2 components to determine mate
/// insert candidacy. Note that we assume base filtration
/// (PCR dups, etc.) has already occurred.
///
bool isMateInsertionEvidenceCandidate(const bam_record& bamRead, const unsigned minMapq);

bool isMateInsertionEvidenceCandidate2(
    const bam_record& bamRead, const bool isSearchForLeftOpen, const bool isSearchForRightOpen);

/// information recorded for reads where we need to grab the mate from a remote locus
///
/// typically these are chimeras with a MAPQ0 mate used to assemble a large insertion
///
struct RemoteReadInfo {
  explicit RemoteReadInfo(const bam_record& bamRead)
    : qname(bamRead.qname()),
      readNo(bamRead.read_no() == 1 ? 2 : 1),
      tid(bamRead.mate_target_id()),
      pos(bamRead.mate_pos() - 1),
      localPos(bamRead.pos() - 1),
      readSize(bamRead.read_size()),
      isLocalFwd(bamRead.is_fwd_strand()),
      isFound(false),
      isUsed(false)
  {
  }

  bool operator<(const RemoteReadInfo& rhs) const
  {
    if (tid < rhs.tid) return true;
    if (tid == rhs.tid) {
      return (pos < rhs.pos);
    }
    return false;
  }

  std::string qname;
  int         readNo;  // this is read number of the target
  int         tid;
  int         pos;
  int         localPos;
  int         readSize;
  bool        isLocalFwd;
  bool        isFound;
  bool        isUsed;
};
