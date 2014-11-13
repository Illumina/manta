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

#include "bed_streamer.hh"
#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "blt_util/seq_util.hh"
#include <cassert>
#include <cstdlib>
#include <sys/stat.h>

#include <iostream>
#include <set>
#include <string>



bool
bed_streamer::
next()
{
    if (_is_stream_end || (nullptr==_hfp) || (nullptr==_titr)) return false;

    char*& record_string(_kstr.s);
    while (true)
    {
        if (tbx_itr_next(_hfp, _tidx, _titr, &_kstr) < 0)
        {
            _is_stream_end=true;
        }
        else
        {
            _is_stream_end=(nullptr == record_string);
        }
        _is_record_set=(! _is_stream_end);
        if (! _is_record_set) break;

        // filter out header for whole file access case:
        if (record_string[0] == '#') continue;

        _record_no++;

        if (! _bedrec.set(record_string))
        {
            log_os << "ERROR: Can't parse BED record: '" << record_string << "'\n";
            exit(EXIT_FAILURE);
        }
        if (! _bedrec.is_valid()) continue;
        break;
    }

    return _is_record_set;
}



void
bed_streamer::
report_state(std::ostream& os) const
{
    const bed_record* bedp(get_record_ptr());

    os << "\tvcf_stream_label: " << name() << "\n";
    if (nullptr != bedp)
    {
        os << "\tbed_stream_record_no: " << record_no() << "\n"
           << "\tbed_record: " << *(bedp) << "\n";
    }
    else
    {
        os << "\tno bed record currently set\n";
    }
}
