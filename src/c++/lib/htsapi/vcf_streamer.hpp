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

#include "hts_streamer.hpp"
#include "vcf_record.hpp"

#include <cassert>

/// \brief Stream records from VCF files.
struct vcf_streamer : public hts_streamer {
  /// \param[in] requireNormalized if true an exception is thrown for any input variant records which are not
  /// normalized (see below for definition)
  ///
  /// A VCF record is considered normalized if it is left-aligned, has a non-zero ref and alt length, and
  /// is parsimonious except for left-side reference padding required to fulfill the non-zero length rule.
  ///
  vcf_streamer(const char* filename, const char* region, const bool requireNormalized = true);

  ~vcf_streamer();

  /// advance to next (normalized) vcf record
  ///
  bool next();

  const vcf_record* get_record_ptr() const
  {
    if (_is_record_set)
      return &_vcfrec;
    else
      return nullptr;
  }

  void report_state(std::ostream& os) const;

  /// provide a BAM header to validate vcf chromosome names against
  ///
  void validateBamHeaderChromSync(const bam_hdr_t& header) const;

  unsigned getSampleCount() const
  {
    assert(nullptr != _hdr);
    return _sampleCount;
  }

  const char* getSampleName(const unsigned sampleIndex) const
  {
    assert(nullptr != _hdr);
    assert(sampleIndex < _sampleCount);
    return _hdr->samples[sampleIndex];
  }

private:
  bcf_hdr_t* _hdr;
  unsigned   _sampleCount;
  vcf_record _vcfrec;
  bool       _requireNormalized;
};
