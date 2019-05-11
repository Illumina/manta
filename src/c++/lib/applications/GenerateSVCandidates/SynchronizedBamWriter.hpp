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

#include <mutex>

#include "boost/noncopyable.hpp"

#include "htsapi/bam_dumper.hpp"

/// Extends the standard BAM writer to allow synchronized writing to the same file from multiple threads
class SynchronizedBamWriter : private boost::noncopyable {
public:
  SynchronizedBamWriter(const char* filename, const bam_hdr_t& header) : m_bamWriter(filename, header) {}

  /// Add another BAM record to the file. File must not be closed.
  void put_record(const bam1_t* brec)
  {
    std::lock_guard<std::mutex> lock(m_writeMutex);
    m_bamWriter.put_record(brec);
  }

private:
  bam_dumper m_bamWriter;
  std::mutex m_writeMutex;
};
