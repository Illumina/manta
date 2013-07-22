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
#include "SvLocusAssembler.hh"

#include <string>
#include <utility>
#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>


// compile with this macro to get verbose output:
//#define DEBUG_ASBL

#ifdef DEBUG_ASBL
#include <iosfwd>
extern std::ostream& dbg_os;
#endif

class SvLocusAssembler {


public:
    SvLocusAssembler() {};

    // a shadow/anomalous/semi-aligning read which stores its position in the cluster and remember if
    // it has been used in an assembly.
    struct ShadowRead
    {
        ShadowRead() : seq("EMPTY"), pos(0), used(false) {}
        ShadowRead(std::string s, unsigned p, bool u) : seq(s), pos(p), used(u) {}

        std::string seq; // sequence
        unsigned pos; // position in cluster
        bool used;		// used in assembly?
    };

    // Vector of shadow reads with their index (i.e. position in the cluster)
    //typedef std::vector<std::pair<unsigned,std::string> > ShadowReadVec;
    typedef std::vector< ShadowRead > ShadowReadVec;

    // maps kmers to positions in read
    typedef boost::unordered_map<std::string,unsigned> str_uint_map_t;

    //
    typedef std::vector<AssembledContig> Assembly;

    /**
     * @brief Performs a simple de-novo assembly using a group of reads
     *
    * Assembles a cluster of shadow reads. This function iterates over a range
    * of word lengths until the first succesful assembly.
     *
    * If unused reads remain, the assembly is re-started using this subset.
     */
    void iterateWordLengths(const SmallAssemblerOptions& ,ShadowReadVec& ,Assembly& );

private:

    /**
     * @brief Performs a simple de-novo assembly using a group of reads
     *
     * We build a hash of the k-mers in the shadow reads. The most
     * frequent k-mer is used as seed and is then gradually extended.
     *
     */
    bool buildContigs(const SmallAssemblerOptions& asmOptions,
                      ShadowReadVec& shadows,
                      const unsigned wordLength,
                      std::vector<AssembledContig>& contigs,
                      unsigned& unused_reads);

    /**
    * Adds base @p base to the end (mode=0) or start (mode=1) of the contig.
     *
     * @return The extended contig.
     */
    std::string addBase(const std::string& contig, const char base, const unsigned int mode);

    /**
    * Returns suffix (mode=0) or prefix (mode=1) of @p contig with length @p length.
    *
    *  @return The suffix or prefix.
    */
    std::string getEnd(const std::string& contig, const unsigned length, const unsigned mode);

    /**
     * Extends the seed contig (aka most frequent k-mer)
    *
    *  @return The extended contig.
    */
    void walk(const SmallAssemblerOptions& asmOptions,
              const std::string& seed,
              const unsigned wordLength,
              const str_uint_map_t& wordHash,
              unsigned& stepsBackward,
              std::string& contig);

};

#endif 
