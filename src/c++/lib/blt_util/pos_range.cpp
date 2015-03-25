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

#include "blt_util/pos_range.hh"

#include <iostream>


// output is always 1-indexed inclusive interval:
//
std::ostream& operator<<(std::ostream& os, const pos_range& pr)
{
    os << "[";
    if (pr.is_begin_pos)
    {
        os << pr.begin_pos+1;
    }
    else
    {
        os << "-inf";
    }
    os << " .. ";
    if (pr.is_end_pos)
    {
        os << pr.end_pos;
    }
    else
    {
        os << "inf";
    }
    os << "]";

    return os;
}
