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

#include "bed_record.hh"
#include "blt_util/parse_util.hh"

#include <cassert>

#include <iostream>



bool
bed_record::
set(const char* s)
{
    static const char sep('\t');
    static const unsigned maxword(3);

    clear();

    line = s;

    // simple tab parse:
    const char* start(s);
    const char* p(start);

    unsigned wordindex(0);
    while (wordindex<maxword)
    {
        if ((*p == sep) || (*p == '\n') || (*p == '\0'))
        {
            switch (wordindex)
            {
            case 0:
                chrom=std::string(start,p-start);
                break;
            case 1:
                begin=illumina::blt_util::parse_int(start);
                assert(start==p);
                break;
            case 2:
                end=illumina::blt_util::parse_int(start);
                assert(start==p);
                break;
            default:
                assert(0);
                break;
            }
            start=p+1;
            wordindex++;
        }
        if ((*p == '\n') || (*p == '\0')) break;
        ++p;
    }

    return (wordindex >= maxword);
}



std::ostream&
operator<<(
    std::ostream& os,
    const bed_record& bedr)
{
    os << bedr.chrom << '\t'
       << bedr.begin << '\t'
       << bedr.end << '\n';

    return os;
}

