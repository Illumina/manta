// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
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

#include "blt_util/known_pos_range2.hh"

#include <iostream>


// output is always 1-indexed inclusive interval:
//
std::ostream& operator<<(std::ostream& os, const known_pos_range2& pr)
{
    os << '['
       << pr.begin_pos()
       << ','
       << pr.end_pos()
       << ')';

    return os;
}
