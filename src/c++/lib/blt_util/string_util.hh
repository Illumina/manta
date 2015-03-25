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

#include <string>
#include <vector>


void
split_string(
    const char* str,
    const char delimiter,
    std::vector<std::string>& v);

/// insert nulls into str to create a vector of c-strs without new allocation
void
destructive_split_string(
    char* str,
    const char delimiter,
    std::vector<const char*>& v);

void
split_string(
    const std::string& str,
    const char delimiter,
    std::vector<std::string>& v,
    const bool isSkipEmpty = false);

/// check for exact match to pattern after delimiting str by delimiter
bool
split_match(
    const std::string& str,
    const char delimiter,
    const char* needle);

