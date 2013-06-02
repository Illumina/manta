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

#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include "blt_util/parse_util.hh"


int
main() {
    using namespace illumina::blt_util;

    static const char two[] = "2";
    parse_int_str(std::string(two));
}
