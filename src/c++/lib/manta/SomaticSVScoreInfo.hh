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

///
/// \author Chris Saunders
///

#pragma once

#include <iosfwd>
#include <string>
#include <vector>


/// consolidate all somatic scoring results applied to an SV candidate
struct SomaticSVScoreInfo
{
    SomaticSVScoreInfo()
    {}

    std::vector<std::string> filters;
};


std::ostream&
operator<<(std::ostream& os, const SomaticSVScoreInfo& ssi);


