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

/// \author Chris Saunders
///

#include "blt_util/blt_exception.hh"

#include <cstdlib>

#include <fstream>
#include <iostream>
#include <sstream>



void
open_ifstream(std::ifstream& ifs,
              const char* filename)
{

    ifs.open(filename);
    if (! ifs)
    {
        std::ostringstream oss;
        oss << "ERROR: Can't open file: " << filename << "\n";
        throw blt_exception(oss.str().c_str());
    }
}
