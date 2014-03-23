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
/// \author Ole Schulz-Trieglaff
///

#pragma once

#include "GSCOptions.hh"

#include "alignment/GlobalAligner.hh"
#include "alignment/GlobalJumpAligner.hh"
#include "alignment/GlobalJumpIntronAligner.hh"
#include "blt_util/bam_header_info.hh"
#include "manta/SVCandidateAssemblyData.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVLocusAssembler.hh"
#include "options/SmallAssemblerOptions.hh"


/// \brief methods to improve low-resolution SVCandidates via assembly and contig alignment
///
struct SVCandidateAssemblyRefiner
{
    SVCandidateAssemblyRefiner(
        const GSCOptions& opt,
        const bam_header_info& header);

    /// \brief add assembly and assembly post-processing data to SV candidate
    ///
    /// \param isRNA[in] if true add intron logic to the contig jump aligner
    /// \param isRNA[in] if true search for insertions which can't be completely assembled
    void
    getCandidateAssemblyData(
        const SVCandidate& sv,
        const SVCandidateSetData& svData,
        const bool isRNA,
        const bool isFindLargeInsertions,
        SVCandidateAssemblyData& assemblyData) const;

private:

    /// large SV assembler
    void
    getJumpAssembly(
        const SVCandidate& sv,
        const bool isRNA,
        const bool isFindLargeInsertions,
        SVCandidateAssemblyData& assemblyData) const;

    /// small SV/indel assembler
    ///
    /// \params[in] isFindLargeInsertions search for insertions which are too large to completely assemble the insert sequence
    ///
    void
    getSmallSVAssembly(
        const SVCandidate& sv,
        const bool isFindLargeInsertions,
        SVCandidateAssemblyData& assemblyData) const;

    //////////////////////////////// data:
    const GSCOptions& _opt;
    const bam_header_info& _header;
    const SVLocusAssembler _smallSVAssembler;
    const SVLocusAssembler _spanningAssembler;
    const GlobalAligner<int> _smallSVAligner;
    const GlobalAligner<int> _largeInsertAligner;
    const GlobalJumpAligner<int> _spanningAligner;
    const GlobalJumpIntronAligner<int> _RNASpanningAligner;
};
