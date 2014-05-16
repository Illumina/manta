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

///
/// \author Bret Barnes, Xiaoyu Chen
///

#pragma once

#include "manta/ReadGroupStatsSet.hh"

#include <string>


void
extractReadGroupStatsFromBam(
    const std::string& statsBamFile,
    ReadGroupStatsSet& rstats,
    const bool isTumor = false);
