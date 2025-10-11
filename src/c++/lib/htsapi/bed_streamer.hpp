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
