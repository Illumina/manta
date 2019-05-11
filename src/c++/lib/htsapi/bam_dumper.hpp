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
