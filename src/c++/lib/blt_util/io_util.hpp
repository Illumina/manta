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

#include <iosfwd>
#include <memory>
#include <mutex>

#include "boost/noncopyable.hpp"

void open_ifstream(std::ifstream& ifs, const char* filename);

/// use this class to set scope specific stream formatting
///
/// see unit test for example usage
///
struct StreamScoper {
  explicit StreamScoper(std::ostream& os);

  ~StreamScoper();

private:
  std::ostream&                  _os;
  std::unique_ptr<std::ofstream> _tmp_os;
};

/// Synchronizes access to a file stream from multiple threads:
class SynchronizedOutputStream : private boost::noncopyable {
public:
  explicit SynchronizedOutputStream(const std::string& outputFile);

  void write(const std::string& msg);

private:
  std::unique_ptr<std::ostream> m_osPtr;
  std::mutex                    m_writeMutex;
};
