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
///

#pragma once

#include "bed_record.hpp"
#include "hts_streamer.hpp"

/// \brief Stream records from BED files.
//
// Example use:
//
// bed_streamer bst("foo.bed.gz","chr20:100-200");
// while (bst.next()) {
//     const bed_record& bre(*(bst.get_record_ptr()));
//     if(bre.end > 100) foo++;
// }
//
struct bed_streamer : public hts_streamer {
  /// \param[in] filename BED filename, must be non-null & non-empty
  ///
  /// \param[in] region Region string in samtools format (eg. "chr2:20-30"), must both be non-null & non-empty
  ///
  /// \param[in] requireNonZeroRegionLength If true, an exception is thrown for any input bed record with
  /// region size of 0 or less, otherwise such records are skipped without error.
  bed_streamer(const char* filename, const char* region, const bool requireNonZeroRegionLength = true);

  /// \brief Advance to next record
  ///
  bool next();

  const bed_record* get_record_ptr() const
  {
    if (_is_record_set)
      return &_bedrec;
    else
      return nullptr;
  }

  void report_state(std::ostream& os) const;

private:
  bed_record _bedrec;
  bool       _requireNonZeroRegionLength;
};
