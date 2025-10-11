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
