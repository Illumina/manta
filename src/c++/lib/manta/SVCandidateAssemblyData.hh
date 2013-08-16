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
/// \author Chris Saunders

#pragma once

#include "assembly/AssembledContig.hh"
#include "alignment/GlobalJumpAligner.hh"

#include <vector>


/// Assembly data pertaining to a specific SV candidate
///
/// Assembly starts from a low-resolution SV candidate. This holds
/// any persistent data related to the assembly process, such as data
/// useful for scoring the candidate.
///
/// As a future design note -- it may be possible that the candidate is broken
/// into multiple candidates during assembly-based hypothesis refinement, so
/// this struct may cover multiple refined candidates (but always only one input
/// candidate.
///
struct SVCandidateAssemblyData
{
    SVCandidateAssemblyData() :
        isBestAlignment(false),
        bestAlignmnetIndex(0)
    {}

    void
    clear()
    {
        contigs.clear();
        alignments.clear();
        isBestAlignment=false;
        bestAlignmnetIndex=0;
    }

    typedef JumpAlignmentResult<int> JumpAlignmentResultType;

    Assembly contigs; ///< assembled contigs for both breakpoints
    std::vector<JumpAlignmentResultType> alignments; ///< contig alignments, one per contig, may be empty

    bool isBestAlignment;        ///< is there a contig/alignment good enough to be used for reporting?
    unsigned bestAlignmnetIndex; ///< if isBestAlignment, which is the best?
};
