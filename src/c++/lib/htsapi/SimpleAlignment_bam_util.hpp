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

///
///
///

#pragma once

#include "blt_util/SimpleAlignment.hpp"
#include "htsapi/bam_record.hpp"

void getAlignment(const bam_record& bamRead, SimpleAlignment& al);

SimpleAlignment getAlignment(const bam_record& bamRead);

/// generate a mate alignment, using the MATE CIGAR (MC) BAM tag if available, else assuming same read length
/// and perfect alignment
SimpleAlignment getKnownOrFakedMateAlignment(const bam_record& bamRead);
