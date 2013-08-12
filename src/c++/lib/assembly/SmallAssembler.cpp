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


#include "assembly/SmallAssembler.hh"

#include "boost/foreach.hpp"
#include "boost/unordered_map.hpp"

#include <cassert>

#include <vector>


// compile with this macro to get verbose output:
//#define DEBUG_ASBL


// stream used by DEBUG_ASBL:
#ifdef DEBUG_ASBL
#include <iostream>
std::ostream& dbg_os(std::cerr);
#endif



// maps kmers to positions in read
typedef boost::unordered_map<std::string,unsigned> str_uint_map_t;
// remembers used reads
typedef boost::unordered_map<std::string,bool> str_bool_map_t;



/**
 * Adds base @p base to the end (mode=0) or start (mode=1) of the contig.
 *
 *	@return The extended contig.
 */
static
std::string
addBase(const std::string& contig,
        const char base,
        const unsigned int mode)
{

    switch (mode)
    {
    case 0:
        return contig + base;
    case 1:
        return base + contig;
    default:
        std::cerr << "ERROR: addBase() : " << contig << " " << base << " " << mode << "\n";
        exit(EXIT_FAILURE);
    }

    // FIXME: Why return empty string?
    return std::string();
}



/**
 * Returns a suffix (mode=0) or prefix (mode=1) of @p contig with length @p length.
 *
 *	@return The suffix or prefix.
 */
static
std::string
getEnd(const std::string& contig,
       const unsigned length,
       const unsigned mode)
{

    const unsigned csize(contig.size());
    assert(length <= csize);

    switch (mode)
    {
    case 0:
        return contig.substr((csize-length),length);
    case 1:
        return contig.substr(0,length);
    default:
        std::cerr << "ERROR:: getEnd() : " << contig << " " << length << " " << mode << "\n";
        exit(EXIT_FAILURE);
    }

    // FIXME: Why return empty string?
    return std::string();
}



/**
 * Extends the seed contig (aka most frequent k-mer)
 *
 */
static
void
walk(const SmallAssemblerOptions& opt,
     const std::string& seed,
     const unsigned wordLength,
     const str_uint_map_t& wordHash,
     std::string& contig)
{
    // we start with the seed
    contig = seed;

    std::set<std::string> seenBefore;	// records k-mers already encountered during extension
    seenBefore.insert(contig);

    const str_uint_map_t::const_iterator whe(wordHash.end());

    // 0 => walk to the right, 1 => walk to the left
    for (unsigned mode(0); mode<2; ++mode)
    {
        while (true)
        {
            const std::string tmp(getEnd(contig,wordLength-1,mode));

#ifdef DEBUG_ASBL
            dbg_os << "# current contig : " << contig << " size : " << contig.size() << "\n"
                   << " getEnd : " << tmp << "\n";
#endif

            if (seenBefore.find(tmp) != seenBefore.end())
            {
#ifdef DEBUG_ASBL
                dbg_os << "Seen word " << tmp << " before on this walk, terminating" << "\n";
#endif
                break;
            }

            seenBefore.insert(tmp);

            unsigned maxOcc(0);
            unsigned totalOcc(0);
            char maxBase('N');

            static const std::string BASES("ACGT");
            BOOST_FOREACH(const char b, BASES)
            {
                const std::string newKey(addBase(tmp,b,mode));
#ifdef DEBUG_ASBL
                dbg_os << "Extending end : base " << b << " " << newKey << "\n";
#endif
                const str_uint_map_t::const_iterator whi(wordHash.find(newKey));
                const unsigned val((whi == whe) ? 0 : whi->second);

                totalOcc += val;
                if (val > maxOcc)
                {
                    maxOcc  = val;
                    maxBase = b;
                }
#ifdef DEBUG_ASBL
                dbg_os << b << " " << newKey << " " << val << " " << "\n";
#endif
            }
#ifdef DEBUG_ASBL
            dbg_os << "Winner is : " << maxBase << " with " << maxOcc << " occurrences." << "\n";
#endif

            if ((maxOcc < opt.minCoverage) ||
                (maxOcc < (1.- opt.maxError)* totalOcc))
            {
#ifdef DEBUG_ASBL
                dbg_os << "Coverage or error rate below threshold.\n"
                       << "maxocc : " << maxOcc << " min coverage: " << _minCoverage << "\n"
                       << "maxocc : " << maxOcc << " error threshold: " << ((1.0-_maxError)* totalOcc) << "\n";
#endif
                break;
            }

#ifdef DEBUG_ASBL
            dbg_os << "Adding base " << contig << " " << maxBase << " " << mode << "\n";
#endif
            contig = addBase(contig,maxBase,mode);
#ifdef DEBUG_ASBL
            dbg_os << "New contig: " << contig << "\n";
#endif
        }
#ifdef DEBUG_ASBL
        dbg_os << "mode change. Current mode " << mode << "\n";
#endif
    }

    return;
}



bool
buildContigs(const SmallAssemblerOptions& opt,
             AssemblyReadMap& reads,
             const unsigned wordLength,
             Assembly& as,
             unsigned& unused_reads)
{
#ifdef DEBUG_ASBL
    dbg_os << "In SVLocusAssembler::buildContig. word length=" << wordLength << "\n";
    for (AssemblyReadMap::const_iterator ct = reads.begin(); ct!=reads.end(); ++ct)
    {
        dbg_os << ct->second.seq << " used=" << ct->second.used << "\n";
    }
#endif

    // a set of read hashes; each read hash stores the starting positions of all kmers in the read
    std::vector<std::pair<std::string,str_uint_map_t> > readHashes;
    // counts the number of occurrences for each kmer in all shadow reads
    str_uint_map_t wordHash;
    // most frequent kmer and its number of occurrences
    unsigned maxOcc = 0;
    std::string maxWord;

    BOOST_FOREACH(const AssemblyReadMap::value_type& contigRead, reads)
    {
        // skip reads used in an previous iteration
        if (contigRead.second.used) continue;

        // stores the index of a kmer in a read sequence
        const std::string& seq(contigRead.second.seq);
        const std::string& qn(contigRead.first);

        const unsigned readLen(seq.size());

        // TODO: (csaunders) should this be reduced from an assertion to something more robust? What if we skipped short reads?
        assert(readLen>=wordLength);

        readHashes.resize(readHashes.size()+1);
        readHashes.back().first  = qn;
        str_uint_map_t& readHash = readHashes.back().second;

        for (unsigned j(0); j<=(readLen-wordLength); ++j)
        {
            const std::string word(seq.substr(j,wordLength));
            if (readHash.find(word) != readHash.end())
            {
                // try again with different k-mer size
#ifdef DEBUG_ASBL
                dbg_os << "word " << word << " repeated in read " << qn << "\n"
                       << "... try again with word length " << wordLength+2
                       << " and " << reads.size() << " reads left.\n";
#endif
                return false;
            }

            // record (0-indexed) start point for word in read
            //cout << "Recording " << word << " at " << j << "\n";
            readHash[word]=j;

            // count occurrences
            ++wordHash[word];
            if (wordHash[word]>maxOcc)
            {
                //cout << "Setting max word to " << maxWord << " " << maxOcc << "\n";
                maxOcc  = wordHash[word];
                maxWord = word;
            } // if (maxOcc)
        }
    }

    if (maxOcc < opt.minCoverage)
    {
#ifdef DEBUG_ASBL
        dbg_os << "Coverage too low : " << maxOcc << " " << _minCoverage << "\n";
#endif
        return false;
    }

#ifdef DEBUG_ASBL
    dbg_os << "Seeding kmer : " << maxWord << "\n";
#endif

    // start initial assembly with most frequent kmer as seed
    AssembledContig ctg;
    walk(opt,maxWord,wordLength,wordHash,ctg.seq);

    // done with this now:
    wordHash.clear();

    const unsigned rh_size(readHashes.size());
    const unsigned contig_size(ctg.seq.size());
    const AssemblyReadMap::const_iterator arm_e = reads.end();

#ifdef DEBUG_ASBL
    dbg_os << "First pass assembly resulted in "
           << ctg.seq << "\n"
           << " with length " << contig_size << ". Input consisted of " << rh_size << " reads.\n";
#endif

    for (unsigned i(0); i<rh_size; ++i) {
        const std::string& qName(readHashes[i].first);
        const str_uint_map_t& readHash(readHashes[i].second);
        const str_uint_map_t::const_iterator rhe(readHash.end());

        // increment number of reads containing the seeding kmer
        if (readHash.find(maxWord) != rhe)
        {
            ++ctg.seedReadCount;
        }

        // store all reads sharing k-mers of the current word length with the contig
        // TODO: check if we still needs this
        for (unsigned j(0); j<=(contig_size-wordLength); ++j)
        {
            const std::string word(ctg.seq.substr(j,wordLength));
            //cout << "Testing word " << word << " " << readNum << "\n";
            //cout << "with counts : " << wordHash[word] << "\n";
            const str_uint_map_t::const_iterator rhi(readHash.find(word));
            if (rhi != rhe)
            {
                AssemblyReadMap::iterator arm_i(reads.find(qName));
                if (arm_i != arm_e)
                {
                    arm_i->second.used = true;
                    --unused_reads;
                }
                else
                {
                    // TODO: is this an error or a warning??
                    std::cerr << "Read " << qName << " not found in hash\n";
                }
            }
        }
    }

    // don't need this anymore:
    readHashes.clear();
#ifdef DEBUG_ASBL
    dbg_os << "final seeding reading count: " << ctg.seedReadCount << "\n";
#endif
    if (ctg.seedReadCount > opt.minSeedReads)
    {
        as.push_back(ctg);
    }
    else
    {
#ifdef DEBUG_ASBL
        dbg_os << "which is below minSeedReadCount of " << opt.minSeedReads << " discarding.\n";
#endif
    }
    return true;
}



void
runSmallAssembler(
    const SmallAssemblerOptions& opt,
    AssemblyReadMap& readMap,
    Assembly& as)
{
#ifdef DEBUG_ASBL
    dbg_os << "SmallAssember: Starting assembly with " << readMap.size() << " readMap.\n";
#endif

    unsigned unusedReadsNow(readMap.size());
    for(unsigned iterations(0); iterations < opt.maxAssemblyIterations; ++iterations)
    {
        const unsigned unusedReadsPrev(unusedReadsNow);
        for (unsigned wl(opt.wordLength); wl<=opt.maxWordLength; wl+=2)
        {
            const bool isAssemblySuccess = buildContigs(opt,readMap,wl,as,unusedReadsNow);
            if (isAssemblySuccess) break;
        }
        //dbg_os << "iter: " << iteration << "unused readMap now: " << unusedReadsNow << " unused readMap previous: " << unusedReadsPrev << "\n";

        if (unusedReadsNow == 0) return;

        // stop if no change in number of unused reads
        if (unusedReadsNow == unusedReadsPrev) return;
    }
#ifdef DEBUG_ASBL
    dbg_os << "SmallAssembler: Reached max number of assembly iterations: " << opt.maxAssemblyIterations << "\n";
#endif
}


