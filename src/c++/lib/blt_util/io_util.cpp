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

#include "blt_util/io_util.hpp"

#include "blt_util/blt_exception.hpp"

#include <cstdlib>

#include <fstream>
#include <iostream>
#include <sstream>

void open_ifstream(std::ifstream& ifs, const char* filename)
{
  ifs.open(filename);
  if (!ifs) {
    std::ostringstream oss;
    oss << "Can't open file: '" << filename << "'";
    throw blt_exception(oss.str().c_str());
  }
}

StreamScoper::StreamScoper(std::ostream& os) : _os(os), _tmp_os(new std::ofstream)
{
  _tmp_os->copyfmt(_os);
}

StreamScoper::~StreamScoper()
{
  _os.copyfmt(*_tmp_os);
}

SynchronizedOutputStream::SynchronizedOutputStream(const std::string& outputFile)
{
  if (outputFile.empty()) {
    std::ostringstream oss;
    oss << "No output file specified to SynchronizedOutputStream";
    throw blt_exception(oss.str().c_str());
  }
  m_osPtr.reset(new std::ofstream(outputFile.c_str()));
  if (!*m_osPtr) {
    std::ostringstream oss;
    oss << "Can't open output file: '" << outputFile << "'";
    throw blt_exception(oss.str().c_str());
  }
}

void SynchronizedOutputStream::write(const std::string& msg)
{
  std::lock_guard<std::mutex> lock(m_writeMutex);
  *m_osPtr << msg;
}
