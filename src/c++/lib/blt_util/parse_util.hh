// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
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

#pragma once

#include "boost/static_assert.hpp"

#include <string>


namespace illumina
{
namespace blt_util
{

/// parse TYPE from char* with several error checks, and advance
/// pointer to end of TYPE input
///
unsigned
parse_unsigned(const char*& s);

int
parse_int(const char*& s);

long
parse_long(const char*& s);

double
parse_double(const char*& s);



/// std::string version of above, no ptr advance obviously. explicit rename
/// of functions guards against unexpected std::string temporaries
///
unsigned
parse_unsigned_str(const std::string& s);

int
parse_int_str(const std::string& s);

long
parse_long_str(const std::string& s);

double
parse_double_str(const std::string& s);



/// template version:
///
template <typename T>
T
parse_type(const char*&)
{
    // no scan_string available for type:
    BOOST_STATIC_ASSERT(sizeof(T)==0);
    return T();
}


template <>
inline
unsigned
parse_type<unsigned>(const char*& s)
{
    return parse_unsigned(s);
}

template <>
inline
int
parse_type<int>(const char*& s)
{
    return parse_int(s);
}

template <>
inline
long
parse_type<long>(const char*& s)
{
    return parse_long(s);
}

template <>
inline
double
parse_type<double>(const char*& s)
{
    return parse_double(s);
}


}
}

