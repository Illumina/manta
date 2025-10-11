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

#include "assembly/AssemblyReadInfo.hpp"
#include "options/IterativeAssemblerOptions.hpp"
#include "options/ReadScannerOptions.hpp"

typedef IterativeAssemblerOptions AssemblerOptions;

/// load all reads form bam into assembly input structure with minimal
/// filtration / input manipulation
///
void extractAssemblyReadsFromBam(
    const ReadScannerOptions& scanOpt,
    const AssemblerOptions&   asmOpt,
    const std::string&        referenceFilename,
    const std::string&        alignmentFilename,
    AssemblyReadInput&        reads);
