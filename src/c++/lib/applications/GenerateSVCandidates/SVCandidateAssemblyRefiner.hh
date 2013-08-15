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
///

# pragma once

#include "GSCOptions.hh"

#include "alignment/GlobalJumpAligner.hh"
#include "assembly/SmallAssemblerOptions.hh"
#include "blt_util/bam_header_info.hh"
#include "manta/SVCandidateAssemblyData.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVLocusAssembler.hh"


/// methods to improve low-resolution SVCandidates via assembly and contig alignment
///
struct SVCandidateAssemblyRefiner
{
    SVCandidateAssemblyRefiner(
        const GSCOptions& opt,
        const bam_header_info& header,
        const SmallAssemblerOptions& assembleOpt,
        const AlignmentScores<int>& alignScore,
        const int jumpScore);

    void
    getCandidateAssemblyData(
        const SVCandidate& sv,
        const SVCandidateSetData& svData,
        SVCandidateAssemblyData& adata) const;

private:
    const GSCOptions& _opt;
    const SVLocusAssembler _svAssembler;
    const bam_header_info& _header;
    const GlobalJumpAligner<int> _aligner;
};
