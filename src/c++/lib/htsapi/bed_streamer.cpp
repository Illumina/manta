// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
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

#include <iostream>
#include <sstream>


bool
bed_streamer::
next()
{
    if (_is_stream_end || (nullptr==_hfp) || (nullptr==_titr)) return false;

    while (true)
    {
        if (tbx_itr_next(_hfp, _tidx, _titr, &_kstr) < 0)
        {
            _is_stream_end=true;
        }
        else
        {
            _is_stream_end=(nullptr == _kstr.s);
        }
        _is_record_set=(! _is_stream_end);
        if (! _is_record_set) break;

        // filter out header for whole file access case:
        if (_kstr.s[0] == '#') continue;

        _record_no++;

        if (! _bedrec.set(_kstr.s))
        {
            std::ostringstream oss;
            oss << "ERROR: Can't parse BED record: '" << _kstr.s << "'\n";
            throw blt_exception(oss.str().c_str());
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

    os << "\tbed_stream_label: " << name() << "\n";
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
