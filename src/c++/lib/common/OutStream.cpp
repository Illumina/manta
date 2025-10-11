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

#include "common/OutStream.hpp"

#include "common/Exceptions.hpp"

#include <fstream>
#include <iostream>

OutStream::OutStream(const std::string& fileName)
  : _isInit(false), _fileName(fileName), _osptr(&std::cout), _ofsptr(new std::ofstream)
{
  if (!_fileName.empty()) {
    std::ofstream test;
    openFile(_fileName, test);
  }
}

// required for unique_ptr:
OutStream::~OutStream() {}

void OutStream::initStream()
{
  if (!_fileName.empty()) {
    openFile(_fileName, *_ofsptr);
    _osptr = _ofsptr.get();
  }
  _isInit = true;
}

void OutStream::openFile(const std::string& filename, std::ofstream& ofs)
{
  ofs.open(filename.c_str());
  if (ofs) return;
  std::ostringstream oss;
  oss << "Can't open output file: '" << filename << "'";
  BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
}
