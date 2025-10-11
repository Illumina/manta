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
/// \author Trevor Ramsay
///

#include "manta/ReadFilter.hpp"

#include "blt_util/align_path.hpp"
#include "common/Exceptions.hpp"
#include "htsapi/align_path_bam_util.hpp"

#include <ostream>

bool isReadFilteredCore(const bam_record& bamRead)
{
  if (bamRead.is_filter())
    return true;
  else if (bamRead.is_dup())
    return true;
  // supplementary reads without SA tag
  else if (bamRead.is_supplementary() && (!bamRead.isSASplit()))
    return true;
  else {
    // hack to work with bwamem '-M' formatting,
    // keep secondary reads when they contain an SA tag
    if (bamRead.is_secondary()) {
      if (!bamRead.isSASplit()) return true;
    }
  }
  return false;
}

bool isReadUnmappedOrFilteredCore(const bam_record& bamRead)
{
  if (isReadFilteredCore(bamRead)) return true;
  return bamRead.is_unmapped();
}
