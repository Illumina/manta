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

// take care of some (mostly C99) functions not available in VS C++
//

#pragma once

#include <string>


#ifdef _MSC_VER
#define snprintf _snprintf
#define strdup _strdup
#endif

#if ((defined(_MSC_VER)) && (_MSC_VER <= 1800))
#undef noexcept
#define noexcept
#endif

double
compat_round(const double x);


const char*
compat_basename(const char* s);


// gets canonical name of paths, but only when these refer to existing items
// returns false on error.
bool
compat_realpath(std::string& path);
