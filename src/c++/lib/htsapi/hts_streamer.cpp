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

#include "hts_streamer.hpp"
#include "common/Exceptions.hpp"

#include <iostream>
#include <sstream>

static const kstring_t kinit = {0, 0, 0};

hts_streamer::hts_streamer(const char* filename, const char* region)
  : _is_record_set(false),
    _is_stream_end(false),
    _record_no(0),
    _stream_name(filename),
    _hfp(nullptr),
    _tidx(nullptr),
    _titr(nullptr),
    _kstr(kinit)
{
  if (!filename) {
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException("hts filename is null ptr"));
  }

  if ('\0' == *filename) {
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException("hts filename is empty string"));
  }

  _hfp = hts_open(filename, "r");
  if (!_hfp) {
    std::ostringstream oss;
    oss << "Failed to open hts file: '" << filename << "'";
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
  }

  _load_index();

  // read only a region of HTS file:
  if (region) {
    resetRegion(region);
  }
}

hts_streamer::~hts_streamer()
{
  if (_titr) tbx_itr_destroy(_titr);
  if (_tidx) tbx_destroy(_tidx);
  if (_hfp) hts_close(_hfp);
  if (_kstr.s) free(_kstr.s);
}

void hts_streamer::resetRegion(const char* region)
{
  if (_titr) tbx_itr_destroy(_titr);

  _titr          = tbx_itr_querys(_tidx, region);
  _is_stream_end = (!_titr);
}

void hts_streamer::_load_index()
{
  if (_tidx) return;

  _tidx = tbx_index_load(name());
  if (!_tidx) {
    std::ostringstream oss;
    oss << "Failed to load index for hts file: '" << name() << "'";
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
  }
}
