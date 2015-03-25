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

#include "compat_util.hh"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <iostream>


#ifdef _WIN32
#include "compat_util_win32_realpath.c"
#endif



bool
compat_realpath(std::string& path)
{
    errno=0;
    const char* newpath(realpath(path.c_str(),NULL));
    if ((NULL==newpath) || (errno!=0))
    {
        if (NULL!=newpath) free((void*)newpath);
        return false;
    }
    path = newpath;
    free((void*)newpath);
    return true;
}



double
compat_round(const double x)
{
    if (x>=0.)
    {
        return std::floor(x+0.5);
    }
    else
    {
        return std::ceil(x-0.5);
    }
}



const char*
compat_basename(const char* str)
{
#ifdef _WIN32
    static const char pathsep('\\');
#else
    static const char pathsep('/');
#endif
    const char* res(strrchr(str,pathsep));
    if (NULL==res) return str;
    return res+1;
}
