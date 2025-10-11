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

#include "sam_util.hpp"

#include <string>

/// \brief Helper for bam file writing
struct bam_dumper {
  bam_dumper(const char* filename, const bam_hdr_t& header);

  /// Dtor closes the file if it is not already closed
  ~bam_dumper() { close(); }

  /// Add another BAM record to the file. File must not be closed.
  void put_record(const bam1_t* brec);

  /// Return name of bam stream
  const char* name() const { return _stream_name.c_str(); }

  /// Close the output bam file
  void close();

private:
  htsFile*         _hfp;
  const bam_hdr_t* _hdr;

  std::string _stream_name;
};
