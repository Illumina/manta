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

#include "EdgeInfo.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVCandidateData.hh"
#include "manta/SVLocusSet.hh"
#include "manta/SomaticSVScoreInfo.hh"

#include <iosfwd>
#include <vector>


void
writeSomaticSVVcfHeader(
    const char* referenceFilename,
    const char* version,
    std::ostream& os);


void
writeSomaticSVToVcf(
    const std::string& referenceFilename,
    const SVLocusSet& set,
    const EdgeInfo& edge,
    const SVCandidateData& svData,
    const unsigned svIndex,
    const SVCandidate& sv,
    const SomaticSVScoreInfo& ssInfo,
    std::ostream& os);
