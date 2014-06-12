// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
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
///

#pragma once


/// Input parameters for SmallAssembler, a simple local de-bruijn graph assembler
///
struct IterativeAssemblerOptions
{
    /// sets reasonable default values for 30x DNA-seq, 100bp reads
    IterativeAssemblerOptions() :
        alphabet("ACGT"),
        minWordLength(41),
        maxWordLength(76),
        wordStepSize(5),
        minContigLength(15),
        //minCoverage(2),
        minCoverage(1),
        minConservativeCoverage(2),
        maxError(0.35),
        minUnusedReads(3),
        minSupportReads(2),
        maxAssemblyCount(10)
    {}

    std::string alphabet; ///< the symbol set used during assembly
    unsigned minWordLength; ///< initial word (kmer) length
    unsigned maxWordLength; ///< max word length
    unsigned wordStepSize;
    unsigned minContigLength; ///< min contig size
    unsigned minCoverage; ///< min. coverage required for contig extension
    unsigned minConservativeCoverage; ///< coverage required for conservative contig sub-range
    double maxError; ///< max error rates allowed during contig extension
    unsigned minUnusedReads; ///< min. number of unused reads to enable search for more contigs
    unsigned minSupportReads; ///< min. number of reads required to start assembly
    unsigned maxAssemblyCount; ///< Max. number of assembly returned for a given set of reads
};
