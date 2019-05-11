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

#include "htsapi/bam_dumper.hpp"

#include "common/Exceptions.hpp"

#include <cassert>

#include <iostream>
#include <sstream>

bam_dumper::bam_dumper(const char* filename, const bam_hdr_t& header) : _hdr(&header), _stream_name(filename)
{
  assert(filename);

  _hfp = hts_open(filename, "wb");

  if (!_hfp) {
    std::ostringstream oss;
    oss << "Failed to open SAM/BAM/CRAM file for writing: '" << filename << "'";
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
  }

  const int headerWriteStatus = sam_hdr_write(_hfp, _hdr);
  if (headerWriteStatus != 0) {
    std::ostringstream oss;
    oss << "Failed to write SAM/BAM/CRAM file header for: '" << filename << "'";
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
  }
}

void bam_dumper::close()
{
  if (!_hfp) return;
  const int closeStatus = hts_close(_hfp);
  if (closeStatus != 0) {
    std::ostringstream oss;
    oss << "Failed to close SAM/BAM/CRAM file. hts_close return code: " << closeStatus << " stream name: '"
        << name() << "'\n";
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
  }
  _hfp = nullptr;
}

void bam_dumper::put_record(const bam1_t* brec)
{
  assert(_hfp);
  const int recordWriteStatus = sam_write1(_hfp, _hdr, brec);
  if (recordWriteStatus < 0) {
    std::ostringstream oss;
    oss << "Failed to write new record to BAM file: '" << name() << "'";
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
  }
}
