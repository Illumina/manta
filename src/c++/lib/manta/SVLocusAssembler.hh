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

#include "applications/GenerateSVCandidates/GSCOptions.hh"
#include "assembly/SmallAssembler.hh"
#include "blt_util/bam_streamer.hh"
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
        const GSCOptions& opt,
        const SmallAssemblerOptions& assmbleOpt);

    /**
     * @brief Performs a de-novo assembly of a set of reads crossing a breakpoint.
     *
     * Iterates over a range of word lengths until the first successful assembly.
     *
     * If unused reads remain, the assembly is re-started using this subset.
     */
    void
    assembleSingleSVBreakend(const SVBreakend& bp,
                             Assembly& as);

    void
    assembleSVBreakends(const SVBreakend& bp1,
                        const SVBreakend& bp2,
                        Assembly& as);

private:
    typedef boost::shared_ptr<bam_streamer> streamPtr;

    typedef std::map<std::string,unsigned> ReadIndexType;

    // Collects the reads crossing an SV breakpoint
    void
    getBreakendReads(
        const SVBreakend& bp,
        ReadIndexType& readIndex,
        AssemblyReadInput& reads);

    const ReadScannerOptions _scanOpt;
    const SmallAssemblerOptions _assembleOpt;

    // functions to detect anomalous read
    SVLocusScanner _readScanner;
    std::vector<streamPtr> _bamStreams;
};

