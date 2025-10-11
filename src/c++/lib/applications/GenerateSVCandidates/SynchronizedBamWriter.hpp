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
