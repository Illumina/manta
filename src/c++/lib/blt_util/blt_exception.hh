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

#pragma once

#include <exception>
#include <string>

/// \brief a minimal exception class
struct blt_exception : public std::exception
{
    blt_exception(const char* s);

    const char* what() const noexcept
    {
        return message.c_str();
    }

    std::string message;
};
