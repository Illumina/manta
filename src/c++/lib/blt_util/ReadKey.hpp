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

#include "blt_util/compat_util.hpp"
#include "htsapi/bam_record.hpp"

#include "boost/utility.hpp"

#include <cstdlib>
#include <cstring>

#include <iosfwd>
#include <string>

/// information required to uniquely identify a read:
///
struct ReadKey {
  ReadKey(const bam_record& br, const bool isCopyPtrs = true)
    : _isCopyPtrs(isCopyPtrs),
      _qname((_isCopyPtrs && (nullptr != br.qname())) ? strdup(br.qname()) : br.qname()),
      _readNo(br.read_no())
  {
    assert(nullptr != _qname);
  }

  ReadKey(const char* initQname, const int initReadNo, const bool isCopyPtrs = true)
    : _isCopyPtrs(isCopyPtrs),
      _qname((_isCopyPtrs && (nullptr != initQname)) ? strdup(initQname) : initQname),
      _readNo(initReadNo)
  {
    assert(nullptr != _qname);
  }

  ReadKey(const ReadKey& rhs)
    : _isCopyPtrs(rhs._isCopyPtrs),
      _qname(_isCopyPtrs ? strdup(rhs._qname) : rhs._qname),
      _readNo(rhs._readNo)
  {
  }

private:
  ReadKey& operator=(const ReadKey& rhs);

public:
  ~ReadKey()
  {
    if (_isCopyPtrs) {
      if (nullptr != _qname) free(const_cast<char*>(_qname));
    }
  }

  int readNo() const { return _readNo; }

  const char* qname() const { return _qname; }

  bool operator<(const ReadKey& rhs) const
  {
    if (readNo() < rhs.readNo()) return true;
    if (readNo() == rhs.readNo()) {
      return (strcmp(qname(), rhs.qname()) < 0);
    }
    return false;
  }

  bool operator==(const ReadKey& rhs) const
  {
    return ((readNo() == rhs.readNo()) && ((0 == strcmp(qname(), rhs.qname()))));
  }

private:
  const bool  _isCopyPtrs;
  const char* _qname;
  const int   _readNo;
};

std::ostream& operator<<(std::ostream& os, const ReadKey& rk);
