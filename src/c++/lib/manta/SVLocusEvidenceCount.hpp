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

#pragma once

/// enumerate evidence type estimated on input for each sample
struct SVLocusEvidenceCount {
  void clear()
  {
    total                     = 0;
    ignored                   = 0;
    anom                      = 0;
    split                     = 0;
    anomAndSplit              = 0;
    indel                     = 0;
    assm                      = 0;
    remoteRecoveryCandidates  = 0;
    splitSupplementarySegment = 0;
  }

  void merge(const SVLocusEvidenceCount& rhs)
  {
    total += rhs.total;
    ignored += rhs.ignored;
    anom += rhs.anom;
    split += rhs.split;
    anomAndSplit += rhs.anomAndSplit;
    indel += rhs.indel;
    assm += rhs.assm;
    remoteRecoveryCandidates += rhs.remoteRecoveryCandidates;
    splitSupplementarySegment += rhs.splitSupplementarySegment;
  }

  template <class Archive>
  void serialize(Archive& ar, const unsigned /* version */)
  {
    ar& total& ignored& anom& split& anomAndSplit& indel& assm& remoteRecoveryCandidates&
                                                                splitSupplementarySegment;
  }

  // using doubles for integral counts here because:
  // (1) counts are potentially very high
  // (2) exact counts don't matter

  /// Total number of non-filtered non-supplement/secondary reads scanned
  double total = 0;

  /// Total number of non-filtered non-supplement/secondary reads ignored for SV purposes
  double ignored = 0;

  /// Total number of non-filtered non-supplement/secondary anomalous reads scanned
  double anom = 0;

  /// Total number of non-filtered non-supplement/secondary split (SA-tag) reads scanned
  double split = 0;

  /// Total number of non-filtered non-supplement/secondary reads scanned which are categorized as BOTH
  /// anomalous and split
  double anomAndSplit = 0;

  /// Total number of non-filtered non-supplement/secondary CIGAR large indel reads scanned
  double indel = 0;

  /// Total number of non-filtered non-supplement/secondary semi-aligned reads scanned
  double assm = 0;

  /// Subset of anom. these are reads which qualify as candidates for remote recovery
  double remoteRecoveryCandidates = 0;

  /// Total number of non-filtered supplement/secondary split read segments
  double splitSupplementarySegment = 0;
};
