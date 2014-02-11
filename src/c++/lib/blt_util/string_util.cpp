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

/// \file

/// \author Chris Saunders
///


#include "string_util.hh"

#include <cstring>

#include <iostream>



void
split_string(const char* str,
             const char delimiter,
             std::vector<std::string>& v)
{

    v.clear();
    while (true)
    {
        const char* next(strchr(str,delimiter));
        if ((NULL == next) || (delimiter == '\0'))
        {
            v.push_back(std::string(str));
            return;
        }
        v.push_back(std::string(str,next-str));
        str = next+1;
    }
}



void
split_string(const std::string& str,
             const char delimiter,
             std::vector<std::string>& v)
{

    v.clear();

    size_t start(0);
    while (true)
    {
        size_t next(str.find(delimiter,start));
        v.push_back(std::string(str.substr(start,next-start)));
        if (next == std::string::npos) return;
        start = next+1;
    }
}



bool
split_match(const std::string& str,
            const char delimiter,
            const char* needle)
{

    size_t start(0);
    while (true)
    {
        size_t next(str.find(delimiter,start));
        if (0 == str.compare(start,next-start,needle)) return true;
        if (next == std::string::npos) break;
        start = next+1;
    }
    return false;
}
