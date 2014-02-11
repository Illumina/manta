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

/// \author Chris Saunders
///

#include "blt_util/seq_printer.hh"

#include <cassert>
#include <cstring>

#include <iostream>



/// pretty print sequence is such a way that it's easy to locate position number
///
void
printSeq(
    const char* seq,
    std::ostream& os)
{
    static const unsigned rowSize(100);
    static const unsigned sectionSize(10);

    assert(NULL != seq);
    const unsigned seqLen(strlen(seq));

    for (unsigned i(0); i<seqLen; ++i)
    {
        if (i)
        {
            if      (0 == (i % rowSize))     os << '\n';
            else if (0 == (i % sectionSize)) os << ' ';
        }
        os << seq[i];
    }
}
