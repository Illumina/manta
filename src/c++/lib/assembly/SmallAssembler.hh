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
// <https://github.com/sequencing/licenses/>
//

///
/// \author Ole Schulz-Trieglaff
///

#pragma once

#include "assembly/AssembledContig.hh"
#include "assembly/AssemblyReadInfo.hh"
#include "options/SmallAssemblerOptions.hh"


/// \brief run a de-bruijn graph assembler intended for small-scale allele discovery
///
/// the assembler iteratively builds multiple contigs through a range of word sizes
///
/// \param[in] opt assembly parameters
/// \param[in] reads the set of reads to use for the assembly
/// \param[out] assembledReadInfo for each read in 'reads', provide information on if and how it was assembled into a contig
/// \param[out] contigs zero to many assembled contigs
///
void
runSmallAssembler(
    const SmallAssemblerOptions& opt,
    const AssemblyReadInput& reads,
    AssemblyReadOutput& assembledReadInfo,
    Assembly& contigs);

