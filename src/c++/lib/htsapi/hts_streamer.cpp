// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "hts_streamer.hh"
#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"

#include <cstdlib>

#include <iostream>



hts_streamer::
hts_streamer(
    const char* filename,
    const char* region) :
    _is_record_set(false),
    _is_stream_end(false),
    _record_no(0),
    _stream_name(filename),
    _hfp(nullptr),
    _tidx(nullptr),
    _titr(nullptr),
    _kstr({0,0,0})
{
    if (nullptr == filename)
    {
        throw blt_exception("hts filename is null ptr");
    }

    if (nullptr == region)
    {
        throw blt_exception("hts region is null ptr");
    }

    if ('\0' == *filename)
    {
        throw blt_exception("hts filename is empty string");
    }

    _hfp = hts_open(filename, "r");
    if (nullptr == _hfp)
    {
        log_os << "ERROR: Failed to open hts file: '" << filename << "'\n";
        exit(EXIT_FAILURE);
    }

    // read from a specific region:
    _tidx = tbx_index_load(filename);
    if (nullptr == _tidx)
    {
        log_os << "ERROR: Failed to load index for hts file: '" << filename << "'\n";
        exit(EXIT_FAILURE);
    }

    // read only a region of HTS file:
    _titr = tbx_itr_querys(_tidx, region);
    if (nullptr == _titr)
    {
        _is_stream_end=true;
    }
}



hts_streamer::
~hts_streamer()
{
    if (nullptr != _titr) tbx_itr_destroy(_titr);
    if (nullptr != _tidx) tbx_destroy(_tidx);
    if (nullptr != _hfp) hts_close(_hfp);
}
