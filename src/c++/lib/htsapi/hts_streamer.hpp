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

#include "bam_util.hpp"
#include "tabix_util.hpp"

#include "boost/utility.hpp"

#include <string>

/// \brief Stream records various htslib file types
struct hts_streamer : private boost::noncopyable {
  /// \param[in] filename (required)
  /// \param[in] region (may be nullptr)
  ///
  /// Note that region must be set before iteration will produce any output
  /// (ie. there is no "whole genome iteration" enabled if no region is provided)
  hts_streamer(const char* filename, const char* region);

  ~hts_streamer();

  const char* name() const { return _stream_name.c_str(); }

  unsigned record_no() const { return _record_no; }

  /// \brief Set new region for indexed file
  ///
  /// \param[in] region htslib-style region string in format: "chromName:beginPos-endPos", cannot be nullptr
  void resetRegion(const char* region);

protected:
  /// \brief Load index if it hasn't been set already
  void _load_index();

  bool        _is_record_set;
  bool        _is_stream_end;
  unsigned    _record_no;
  std::string _stream_name;

  htsFile*   _hfp;
  tbx_t*     _tidx;
  hts_itr_t* _titr;
  kstring_t  _kstr;
};
