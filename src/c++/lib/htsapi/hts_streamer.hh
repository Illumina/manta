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

#pragma once

#include "bam_util.hh"
#include "tabix_util.hh"

#include "boost/utility.hpp"


struct hts_streamer : private boost::noncopyable
{
    hts_streamer(
        const char* filename,
        const char* region);

    ~hts_streamer();

protected:
    bool _is_record_set;
    bool _is_stream_end;
    unsigned _record_no;
    std::string _stream_name;

    htsFile* _hfp;
    tbx_t* _tidx;
    hts_itr_t* _titr;
    kstring_t _kstr;
};
