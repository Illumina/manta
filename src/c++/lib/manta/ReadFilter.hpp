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

/// \file
/// \brief Filtering logic common to multiple modules.
/// \author Trevor Ramsay
///

#pragma once

#include "htsapi/bam_record.hpp"

/// This predicate runs isReadFiltered without the mapq components
/// \param bamRead The read to test.
/// \return True if the read is filtered out based on core alignment flags.
bool isReadFilteredCore(const bam_record& bamRead);

/// Test if the read is unmapped or is true in isReadFilteredCore
/// \param bamRead The read to test.
/// \return True if the read is filtered out based on the core alignment flags, or is unmapped.
bool isReadUnmappedOrFilteredCore(const bam_record& bamRead);
