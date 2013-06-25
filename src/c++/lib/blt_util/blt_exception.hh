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
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file

/// \author Chris Saunders
///
#ifndef __BLT_EXCEPTION_HH
#define __BLT_EXCEPTION_HH

#include <exception>
#include <string>

/// \brief a minimal exception class
struct blt_exception : public std::exception
{

    blt_exception(const char* s);

    ~blt_exception() throw() {}

    const char* what() const throw()
    {
        return message.c_str();
    }

    std::string message;
};


#endif

