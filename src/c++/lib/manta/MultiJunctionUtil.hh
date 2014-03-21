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
/// \author Chris Saunders
///

#pragma once

#include "manta/SVCandidate.hh"
#include "manta/SVMultiJunctionCandidate.hh"

#include <vector>


void
findMultiJunctionCandidates(
    const std::vector<SVCandidate>& svs,
    std::vector<SVMultiJunctionCandidate>& mjSVs);
