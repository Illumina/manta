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
/// \author Ole Schulz-Trieglaff
/// \author Chris Saunders
///

#include "RemoteMateReadUtil.hpp"

#include <cstdlib>

bool isMateInsertionEvidenceCandidate(const bam_record& bamRead, const unsigned minMapq)
{
  if (!bamRead.is_paired()) return false;
  if (bamRead.isNonStrictSupplement()) return false;
  if (bamRead.is_unmapped() || bamRead.is_mate_unmapped()) return false;

  if (bamRead.map_qual() < minMapq) return false;

  if (bamRead.target_id() < 0) return false;
  if (bamRead.mate_target_id() < 0) return false;

  if (bamRead.target_id() != bamRead.mate_target_id()) return true;

  /// TODO: better candidate definition based on fragment size distro:
  static const int minSize(10000);
  return (std::abs(bamRead.pos() - bamRead.mate_pos()) >= minSize);
}

bool isMateInsertionEvidenceCandidate2(
    const bam_record& bamRead, const bool isSearchForLeftOpen, const bool isSearchForRightOpen)
{
  if ((!isSearchForLeftOpen) && (!bamRead.is_fwd_strand())) return false;
  if ((!isSearchForRightOpen) && bamRead.is_fwd_strand()) return false;
  return true;
}
