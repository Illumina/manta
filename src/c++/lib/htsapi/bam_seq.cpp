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

/// \file

/// \author Chris Saunders
///

#include "htsapi/bam_seq.hh"

#include <iostream>


std::ostream&
operator<<(std::ostream& os,
           const bam_seq_base& bs)
{

    const unsigned rs(bs.size());
    for (unsigned i(0); i<rs; ++i)
    {
        os << bs.get_char(i);
    }
    return os;
}
