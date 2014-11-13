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

#include "bed_record.hh"
#include "hts_streamer.hh"


struct bed_streamer : private hts_streamer
{
    bed_streamer(
        const char* filename,
        const char* region) :
        hts_streamer(filename,region)
    {}

    /// advance to next record
    ///
    bool next();

    const bed_record*
    get_record_ptr() const
    {
        if (_is_record_set) return &_bedrec;
        else                return nullptr;
    }

    const char* name() const
    {
        return _stream_name.c_str();
    }

    unsigned record_no() const
    {
        return _record_no;
    }

    void report_state(std::ostream& os) const;

private:
    bed_record _bedrec;
};
