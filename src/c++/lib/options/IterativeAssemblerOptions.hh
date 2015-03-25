// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Xiaoyu Chen
/// \author Ole Schulz-Trieglaff
///

#pragma once


/// Input parameters for IterativeAssembler
///
struct IterativeAssemblerOptions
{
    /// the symbol set used during assembly
    std::string alphabet = "ACGT";

    /// initial word (kmer) length
    unsigned minWordLength = 41;

    unsigned maxWordLength = 76;
    unsigned wordStepSize = 5;
    unsigned minContigLength = 15;

    /// min. coverage required for contig extension
    unsigned minCoverage = 1;

    /// coverage required for conservative contig sub-range
    unsigned minConservativeCoverage = 2;

    /// max error rates allowed during contig extension
    double maxError = 0.35;

    /// min. number of unused reads to enable search for more contigs
    unsigned minUnusedReads = 3;

    /// min. number of reads required to start assembly
    unsigned minSupportReads = 2;

    /// Max. number of assembly returned for a given set of reads
    unsigned maxAssemblyCount = 10;
};
