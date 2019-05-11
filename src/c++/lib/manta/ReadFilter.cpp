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
