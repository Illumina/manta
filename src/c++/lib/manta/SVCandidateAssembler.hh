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

#include "applications/GenerateSVCandidates/GSCOptions.hh"
#include "assembly/IterativeAssembler.hh"
#include "assembly/SmallAssembler.hh"
#include "htsapi/bam_streamer.hh"
#include "manta/ChromDepthFilterUtil.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVCandidateAssemblyData.hh"
#include "manta/SVCandidateSetData.hh"
#include "manta/SVLocusScanner.hh"

#include <vector>


//#define ITERATIVE_ASSEMBLER
#ifdef ITERATIVE_ASSEMBLER
typedef IterativeAssemblerOptions AssemblerOptions;
#else
typedef SmallAssemblerOptions AssemblerOptions;
#endif

/// Assembles SV-candidate reads for single and paired SVBreakend objects
///
struct SVCandidateAssembler
{
    SVCandidateAssembler(
        const ReadScannerOptions& scanOpt,
        const AssemblerOptions& assembleOpt,
        const AlignmentFileOptions& alignFileOpt,
        const std::string& statsFilename,
        const std::string& chromDepthFilename,
        const bam_header_info& bamHeader);

    /**
     * @brief Performs a de-novo assembly of a set of reads crossing a breakpoint.
     *
     * Iterates over a range of word lengths until the first successful assembly.
     *
     * If unused reads remain, the assembly is re-started using this subset.
     */
    void
    assembleSingleSVBreakend(
        const SVBreakend& bp,
        const reference_contig_segment& refSeq,
        const bool isSearchRemoteInsertionReads,
        RemoteReadCache& remoteReads,
        Assembly& as) const;

    void
    assembleSVBreakends(
        const SVBreakend& bp1,
        const SVBreakend& bp2,
        const bool isBp1Reversed,
        const bool isBp2Reversed,
        const reference_contig_segment& refSeq1,
        const reference_contig_segment& refSeq2,
        Assembly& as) const;

    const AssemblerOptions&
    getAssembleOpt() const
    {
        return _assembleOpt;
    }


    typedef std::map<std::string,unsigned> ReadIndexType;

private:
    typedef std::shared_ptr<bam_streamer> streamPtr;

    /// Collects the reads crossing an SV breakpoint and adds them to reads
    ///
    /// \param[in] isReversed if true revcomp all reads on input
    /// \param[in] isSearchRemoteInsertionReads if true search the remote end of chimeric pairs for MAPQ0 insertion support
    /// \param[out] remoteReadsCache stores any discovered remote reads so that these can be reused during scoring
    void
    getBreakendReads(
        const SVBreakend& bp,
        const bool isReversed,
        const reference_contig_segment& refSeq,
        const bool isSearchRemoteInsertionReads,
        RemoteReadCache& remoteReadsCache,
        ReadIndexType& readIndex,
        AssemblyReadInput& reads) const;

    const ReadScannerOptions _scanOpt;
    const AssemblerOptions _assembleOpt;
    const std::vector<bool> _isAlignmentTumor;
    const ChromDepthFilterUtil _dFilter;
    const ChromDepthFilterUtil _dFilterRemoteReads;

    // contains functions to detect/classify anomalous reads
    SVLocusScanner _readScanner;
    std::vector<streamPtr> _bamStreams;
};
