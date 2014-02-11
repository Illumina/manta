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

#pragma once

#include <string>


/// check if input file exists and is usable as
/// input, if so canonicalize the name
///
/// In case of error return true and provide error
/// message
bool
checkStandardizeInputFile(
    std::string& filename,
    const char* fileLabel,
    std::string& errorMsg);
