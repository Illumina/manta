//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
