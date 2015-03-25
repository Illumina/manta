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

#include <iosfwd>
#include <memory>


void
open_ifstream(
    std::ifstream& ifs,
    const char* filename);


/// use this class to set scope specific stream formatting
///
/// see unit test for example usage
///
struct StreamScoper
{
    StreamScoper(
        std::ostream& os);

    ~StreamScoper();

private:
    std::ostream& _os;
    std::unique_ptr<std::ofstream> _tmp_os;
};
