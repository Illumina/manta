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
