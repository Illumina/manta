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

#include "blt_util/seq_util.hh"

#include <iosfwd>
#include <string>
#include <vector>


struct bed_record
{
    bed_record()
    {
        clear();
    }

    // set record from record string s, return false on error
    bool set(const char* s);

    void clear()
    {
        chrom.clear();
        begin=0;
        end=0;
        line=nullptr;
    }

    bool
    is_valid() const
    {
        return (begin <= end);
    }

    std::string chrom;
    int begin = 0;
    int end = 0;
    const char* line = nullptr;
};


std::ostream& operator<<(std::ostream& os, const bed_record& bedr);

