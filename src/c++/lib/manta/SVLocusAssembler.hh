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
#include "assembly/SmallAssembler.hh"
#include "blt_util/bam_streamer.hh"
#include "manta/ChromDepthFilterUtil.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVCandidateSetData.hh"
#include "manta/SVLocusScanner.hh"

#include "boost/shared_ptr.hpp"

#include <vector>


/// Assembles SV-candidate reads for single and paired SVBreakend objects
///
struct SVLocusAssembler
{
    SVLocusAssembler(
        const ReadScannerOptions& scanOpt,
        const SmallAssemblerOptions& assembleOpt,
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
        Assembly& as) const;

    void
    assembleSVBreakends(const SVBreakend& bp1,
                        const SVBreakend& bp2,
                        const bool isBp1Reversed,
                        const bool isBp2Reversed,
                        const reference_contig_segment& refSeq1,
                        const reference_contig_segment& refSeq2,
                        Assembly& as) const;

    const SmallAssemblerOptions&
    getAssembleOpt() const
    {
        return _assembleOpt;
    }

private:
    typedef boost::shared_ptr<bam_streamer> streamPtr;

    typedef std::map<std::string,unsigned> ReadIndexType;

    /// Collects the reads crossing an SV breakpoint and adds them to reads
    ///
    /// \param[in] isReversed if true revcomp all reads on input
    /// \param[in] isSearchRemoteInsertionReads if true search the remote end of chimeric pairs for MAPQ0 insertion support
    void
    getBreakendReads(
        const SVBreakend& bp,
        const bool isReversed,
        const reference_contig_segment& refSeq,
        const bool isSearchRemoteInsertionReads,
        ReadIndexType& readIndex,
        AssemblyReadInput& reads) const;

    const ReadScannerOptions _scanOpt;
    const SmallAssemblerOptions _assembleOpt;
    const std::vector<bool> _isAlignmentTumor;
    const ChromDepthFilterUtil _dFilter;

    // contains functions to detect/classify anomalous reads
    SVLocusScanner _readScanner;
    std::vector<streamPtr> _bamStreams;
};
