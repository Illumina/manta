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
/// \author Xiaoyu Chen
///

#pragma once

#include "assembly/AssembledContig.hpp"
#include "assembly/AssemblyReadInfo.hpp"
#include "options/IterativeAssemblerOptions.hpp"

/// \brief run a de-bruijn graph assembler intended for small-scale allele discovery
///
/// the assembler iteratively builds multiple contigs through a range of word sizes
///
/// \param[in] opt assembly parameters
///
/// \param[in] reads the set of reads to use for the assembly
///
/// \param[out] assembledReadInfo for each read in 'reads', provide information on if and how it was assembled
/// into a contig
///
/// \param[out] contigs zero to many assembled contigs
///
void runIterativeAssembler(
    const IterativeAssemblerOptions& opt,
    AssemblyReadInput&               reads,
    AssemblyReadOutput&              assembledReadInfo,
    Assembly&                        contigs);
