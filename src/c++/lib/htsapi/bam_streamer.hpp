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

#include "htsapi/bam_record.hpp"
#include "htsapi/sam_util.hpp"

#include "boost/utility.hpp"

#include <iosfwd>
#include <string>

/// Interface for any object which provides current record and file position for error reporting purposes
struct stream_state_reporter {
  virtual void report_state(std::ostream& /*os*/) const {}

  virtual ~stream_state_reporter();
};

/// Stream bam records from CRAM/BAM/SAM files. For CRAM/BAM
/// files you can run an indexed stream from a specific genome region.
///
//
// Example use:
// bam_streamer stream("sample1.bam","hg19.fasta");
// stream.resetRegion("chr1:1000000-2000000");
// while (stream.next()) {
//     const bam_record& read(*(stream.get_record_ptr()));
//     if(read.is_unmapped()) unmappedCount++;
// }
//
struct bam_streamer : public stream_state_reporter, public boost::noncopyable {
  /// \param filename CRAM/BAM input file
  ///
  /// \param referenceFilename Corresponding reference file. nullptr can be given here to indicate that the
  /// reference is not being provided, but many CRAM files cannot be read in this case.
  ///
  /// \param region Restrict the stream to iterate through a specific region. The BAM/CRAM input file must be
  /// indexed for this option to work. If 'region' is not provided, the stream is configured to iterate
  /// through the entire alignment file.
  bam_streamer(const char* filename, const char* referenceFilename, const char* region = nullptr);

  ~bam_streamer() override;

  /// \brief Set new region to iterate over, this will fail if the alignment file is not indexed
  ///
  /// \param region htslib-style region string in format: "chromName:beginPos-endPos", cannot be nullptr
  void resetRegion(const char* region);

  /// \brief Set new region to iterate over, this will fail if the alignment file is not indexed
  ///
  /// \param referenceContigId htslib zero-indexed contig id
  /// \param beginPos start position (zero-indexed, closed)
  /// \param endPos end position (zero-indexed, closed)
  void resetRegion(int referenceContigId, int beginPos, int endPos);

  bool next();

  const bam_record* get_record_ptr() const
  {
    if (_is_record_set)
      return &_brec;
    else
      return nullptr;
  }

  const char* name() const { return _stream_name.c_str(); }

  unsigned record_no() const { return _record_no; }

  void report_state(std::ostream& os) const override;

  const char* target_id_to_name(const int32_t tid) const;

  int32_t target_name_to_id(const char* seq_name) const;

  const bam_hdr_t& get_header() const { return *(_hdr); }

private:
  void _load_index();

  bool       _is_record_set;
  htsFile*   _hfp;
  bam_hdr_t* _hdr;
  hts_idx_t* _hidx;
  hts_itr_t* _hitr;
  bam_record _brec;

  // track for debug only:
  unsigned    _record_no;
  std::string _stream_name;
  bool        _is_region;
  std::string _region;
};
