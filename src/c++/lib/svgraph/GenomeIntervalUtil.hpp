//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
//
//

#pragma once

#include "svgraph/GenomeInterval.hpp"

#include <string>
#include <vector>

struct bam_header_info;

/// given a collection of genome intervals, reduce down to the minimum non-overlapping set:
///
/// \returns a vector with size equal to the input vector, containing a mapping of the input interval index to
/// the output interval index
///
std::vector<unsigned> intervalCompressor(std::vector<GenomeInterval>& intervals);

/// \brief Build a new GenomeInterval from a samtools-style region string
GenomeInterval convertSamtoolsRegionToGenomeInterval(
    const bam_header_info& bamHeader, const std::string& region);
