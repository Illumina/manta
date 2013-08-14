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



/**
 * Adds base @p base to the end (isEnd is true) or start (otherwise) of the contig.
 *
 *	@return The extended contig.
 */
static
std::string
addBase(const std::string& contig,
        const char base,
        const bool isEnd)
{
    if(isEnd) return contig + base;
    else      return base + contig;
}



/**
 * Returns a suffix (isEnd is true) or prefix (otherwise) of @p contig with length @p length.
 *
 *	@return The suffix or prefix.
 */
static
std::string
getEnd(const std::string& contig,
       const unsigned length,
       const bool isEnd)
{

    const unsigned csize(contig.size());
    assert(length <= csize);

    if(isEnd) return contig.substr((csize-length),length);
    else      return contig.substr(0,length);
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
     const str_uint_map_t& wordCount,
     std::string& contig)
{
    // we start with the seed
    contig = seed;

    std::set<std::string> seenBefore;	// records k-mers already encountered during extension
    seenBefore.insert(contig);

    const str_uint_map_t::const_iterator whe(wordCount.end());

    // 0 => walk to the right, 1 => walk to the left
    for (unsigned mode(0); mode<2; ++mode)
    {
        const bool isEnd(mode==0);

        while (true)
        {
            const std::string tmp(getEnd(contig, wordLength-1, isEnd));

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

            unsigned maxBaseCount(0);
            unsigned totalBaseCount(0);
            char maxBase('N');

            static const std::string BASES("ACGT");
            BOOST_FOREACH(const char b, BASES)
            {
                const std::string newKey(addBase(tmp,b,isEnd));
#ifdef DEBUG_ASBL
                dbg_os << "Extending end : base " << b << " " << newKey << "\n";
#endif
                const str_uint_map_t::const_iterator whi(wordCount.find(newKey));
                const unsigned val((whi == whe) ? 0 : whi->second);

                totalBaseCount += val;
                if (val > maxBaseCount)
                {
                    maxBaseCount  = val;
                    maxBase = b;
                }
            }
#ifdef DEBUG_ASBL
            dbg_os << "Winner is : " << maxBase << " with " << maxBaseCount << " occurrences." << "\n";
#endif

            if ((maxBaseCount < opt.minCoverage) ||
                (maxBaseCount < (1.- opt.maxError)* totalBaseCount))
            {
#ifdef DEBUG_ASBL
                dbg_os << "Coverage or error rate below threshold.\n"
                       << "maxBaseCount : " << maxBaseCount << " minCverage: " << opt.minCoverage << "\n"
                       << "error threshold: " << ((1.0-opt.maxError)* totalBaseCount) << "\n";
#endif
                break;
            }

#ifdef DEBUG_ASBL
            dbg_os << "Adding base " << contig << " " << maxBase << " " << mode << "\n";
#endif
            contig = addBase(contig,maxBase,isEnd);
#ifdef DEBUG_ASBL
            dbg_os << "New contig: " << contig << "\n";
#endif
        }
#ifdef DEBUG_ASBL
        dbg_os << "mode change. Current mode " << mode << "\n";
#endif
    }
}



bool
buildContigs(
    const SmallAssemblerOptions& opt,
    const AssemblyReadInput& reads,
    AssemblyReadOutput& readInfo,
    const unsigned wordLength,
    Assembly& as,
    unsigned& unusedReads)
{
    const unsigned readCount(reads.size());

#ifdef DEBUG_ASBL
    dbg_os << "In SVLocusAssembler::buildContig. word length=" << wordLength << "\n";
    for (unsigned readIndex(0); readIndex<readCount; ++readIndex)
    {
        dbg_os << reads[readIndex] << " used=" << readInfo[readIndex].isUsed << "\n";
    }
#endif

    // a set of read hashes; each read hash stores the starting positions of all kmers in the read
    std::vector<str_uint_map_t> readWordOffsets(readCount);

    // counts the number of occurrences for each kmer in all reads
    str_uint_map_t wordCount;

    // most frequent kmer and its number of occurrences
    unsigned maxWordCount(0);
    std::string maxWord;

    for (unsigned readIndex(0); readIndex<readCount; ++readIndex)
    {
        const AssemblyReadInfo& rinfo(readInfo[readIndex]);

        // skip reads used in a previous iteration
        if (rinfo.isUsed) continue;

        // stores the index of a kmer in a read sequence
        const std::string& seq(reads[readIndex]);
        const unsigned readLen(seq.size());

        // TODO: (csaunders) should this be reduced from an assertion to something more robust? What if we skipped short reads instead of asserting?
        assert(readLen>=wordLength);

        str_uint_map_t& readWordOffset(readWordOffsets[readIndex]);

        for (unsigned j(0); j<=(readLen-wordLength); ++j)
        {
            const std::string word(seq.substr(j,wordLength));
            if (readWordOffset.find(word) != readWordOffset.end())
            {
                // try again with different k-mer size
#ifdef DEBUG_ASBL
                dbg_os << "word " << word << " repeated in read " << readIndex << "\n";
#endif
                return false;
            }

            // record (0-indexed) start point for word in read
            //cout << "Recording " << word << " at " << j << "\n";
            readWordOffset[word]=j;

            // count occurrences
            ++wordCount[word];
            if (wordCount[word]>maxWordCount)
            {
                //cout << "Setting max word to " << maxWord << " " << maxOcc << "\n";
                maxWordCount  = wordCount[word];
                maxWord = word;
            }
        }
    }

    if (maxWordCount < opt.minCoverage)
    {
#ifdef DEBUG_ASBL
        dbg_os << "Coverage too low : " << maxWordCount << " " << opt.minCoverage << "\n";
#endif
        return false;
    }

#ifdef DEBUG_ASBL
    dbg_os << "Seeding kmer : " << maxWord << "\n";
#endif

    // start initial assembly with most frequent kmer as seed
    AssembledContig contig;
    walk(opt,maxWord,wordLength,wordCount,contig.seq);

    // done with this now:
    wordCount.clear();

    const unsigned contigSize(contig.seq.size());

#ifdef DEBUG_ASBL
    dbg_os << "First pass assembly resulted in "
           << contig.seq << "\n"
           << " with length " << contigSize << ". Input consisted of " << readCount << " reads.\n";
#endif

    // increment number of reads containing the seeding kmer
    for (unsigned readIndex(0); readIndex<readCount; ++readIndex)
    {
        const str_uint_map_t& readWordOffset(readWordOffsets[readIndex]);
        if (readWordOffset.count(maxWord)) ++contig.seedReadCount;
    }

#ifdef DEBUG_ASBL
    dbg_os << "final seeding reading count: " << contig.seedReadCount << "\n";
#endif
    if (contig.seedReadCount < opt.minSeedReads)
    {
#ifdef DEBUG_ASBL
        dbg_os << "which is below minSeedReadCount of " << opt.minSeedReads << " discarding.\n";
#endif
        return false;
    }

    // finally -- set isUsed and decrement unusedReads
    for (unsigned readIndex(0); readIndex<readCount; ++readIndex)
    {
        const str_uint_map_t& readWordOffset(readWordOffsets[readIndex]);
        AssemblyReadInfo& rinfo(readInfo[readIndex]);

        if(rinfo.isUsed) continue;

        // store all reads sharing k-mers of the current word length with the contig
        // TODO: check if we still needs this
        for (unsigned j(0); j<=(contigSize-wordLength); ++j)
        {
            const std::string word(contig.seq.substr(j,wordLength));
            //cout << "Testing word " << word << " " << readNum << "\n";
            //cout << "with counts : " << wordCount[word] << "\n";
            if(readWordOffset.count(word))
            {
                rinfo.isUsed = true;
                rinfo.contigId = as.size();
                --unusedReads;
                break;
            }
        }
    }

    // don't need this anymore:
    readWordOffsets.clear();

    as.push_back(contig);
    return true;
}



void
runSmallAssembler(
    const SmallAssemblerOptions& opt,
    const AssemblyReadInput& reads,
    AssemblyReadOutput& assembledReadInfo,
    Assembly& as)
{
#ifdef DEBUG_ASBL
    dbg_os << "SmallAssember: Starting assembly with " << reads.size() << " read.\n";
#endif

    assembledReadInfo.clear();
    as.clear();

    assembledReadInfo.resize(reads.size());

    unsigned unusedReadsNow(reads.size());
    for(unsigned iterations(0); iterations < opt.maxAssemblyIterations; ++iterations)
    {
        const unsigned unusedReadsPrev(unusedReadsNow);
        for (unsigned wordLength(opt.minWordLength); wordLength<=opt.maxWordLength; wordLength+=2)
        {
            const bool isAssemblySuccess = buildContigs(opt, reads, assembledReadInfo, wordLength, as, unusedReadsNow);
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

