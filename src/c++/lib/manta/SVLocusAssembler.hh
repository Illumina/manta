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

#include "AssembledContig.hh"
#include "SVCandidate.hh"
#include "SVCandidateData.hh"
#include "SVLocusScanner.hh"

#include "blt_util/bam_streamer.hh"
//#include "blt_util/log.hh"
#include "blt_util/align_path.hh"

#include "applications/GenerateSVCandidates/GSCOptions.hh"

#include <string>
#include <utility>
#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>

// TODO: change hard-coded settings (Options class ?, check with Chris)


// compile with this macro to get verbose output:
//#define DEBUG_ASBL

#ifdef DEBUG_ASBL
#include <iosfwd>
extern std::ostream& dbg_os;
#endif

class SVLocusAssembler
{

public:

    struct AssemblyRead
    {
        AssemblyRead()
            : seq(""),used(false) {}

        AssemblyRead(std::string s,
                     bool u)
            : seq(s), used(u) {}

        std::string seq;
        bool used;
    };

    // maps kmers to positions in read
    typedef boost::unordered_map<std::string,unsigned> str_uint_map_t;
    // remembers used reads
    typedef boost::unordered_map<std::string,bool> str_bool_map_t;

    //typedef std::vector<AssemblyRead> AssemblyReadMap;
    typedef std::map<std::string,AssemblyRead> AssemblyReadMap;

    typedef boost::shared_ptr<bam_streamer> streamPtr;


    SVLocusAssembler(const GSCOptions& opt);

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

    void
    iterateAssembly(AssemblyReadMap& map,
                    Assembly& as);

    // Collects the reads crossing an SV breakpoint
    void
    getBreakendReads(const SVBreakend& bp,
                     AssemblyReadMap& reads);

    /**
     * @brief Performs a de-novo assembly of a read group.
     *
     * Build a hash of the k-mers in the shadow reads. The most
     * frequent k-mer is used as seed and extended in a greedy fashion.
     */
    bool
    buildContigs(AssemblyReadMap& reads,
                 const unsigned wordLength,
                 Assembly& contigs,
                 unsigned& unused_reads);

    /**
    * Adds base @p base to the end (mode=0) or start (mode=1) of the contig.
     *
     * @return The extended contig.
     */
    std::string
    addBase(const std::string& contig,
            const char base,
            const unsigned int mode);

    /**
    * Returns suffix (mode=0) or prefix (mode=1) of @p contig with length @p length.
    *
    *  @return The suffix or prefix.
    */
    std::string
    getEnd(const std::string& contig,
           const unsigned length,
           const unsigned mode);

    /**
     * Extends the seed contig (aka most frequent k-mer)
    *
    *  @return The extended contig.
    */
    void
    walk(const std::string& seed,
         const unsigned wordLength,
         const str_uint_map_t& wordHash,
         std::string& contig);

    const ReadScannerOptions _scanOpt;
    // functions to detect anomalous read
    SVLocusScanner _readScanner;
    std::vector<streamPtr> _bamStreams;

    //  initial word (kmer) length
    unsigned _wordLength;
    // max word length
    unsigned _maxWordLength;
    // min contig size
    unsigned _minContigLength;
    // min. coverage required for contig extension
    unsigned _minCoverage;
    // max error rates allowed during contig extension
    double  _maxError;
    // min. number of reads required to start assembly
    unsigned _minSeedReads;
    // Max. number of assembly iterations for a cluster before we give up
    unsigned _maxAssemblyIterations;
};

