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
// <https://github.com/downloads/sequencing/licenses/>.
//

///
/// \author Ole Schulz-Trieglaff
///

#pragma once

struct SmallAssemblerOptions
{

    SmallAssemblerOptions() :
        // reasonable default values for 30x and 100bp reads
        wordLength(37), maxWordLength(65), minContigLength(15),
        minCoverage(1), maxError(0.2), minSeedReads(2),
        maxAssemblyIterations(50), contigRangeWobble(200) {}

    //  initial word (kmer) length
    unsigned wordLength;
    // max word length
    unsigned maxWordLength;
    // min contig size
    unsigned minContigLength;
    // min. coverage required for contig extension
    unsigned minCoverage;
    // max error rates allowed during contig extension
    double  maxError;
    // min. number of reads required to start assembly
    unsigned minSeedReads;
    // Max. number of assembly iterations for a cluster before we give up
    unsigned maxAssemblyIterations;
    // Is added to either side of contig alignment range
    unsigned contigRangeWobble;
};

