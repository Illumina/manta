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
