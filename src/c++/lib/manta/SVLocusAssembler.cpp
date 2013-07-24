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


#include "SVLocusAssembler.hh"

#include <cassert>
#include <cstdlib>

#include <iostream>
#include <set>
#include <string>

using namespace std;

// stream used by DEBUG_ASBL:
#ifdef DEBUG_ASBL
ostream& dbg_os(cerr);
#endif

/**
 * Adds base @p base to the end (mode=0) or start (mode=1) of the contig.
 *
 *	@return The extended contig.
 */
string SVLocusAssembler::addBase(const string& contig,
                                 const char base,
                                 const unsigned int mode) {

    switch (mode) {
    case 0:
        return contig + base;
    case 1:
        return base + contig;
    default:
        cerr << "ERROR: addBase() : " << contig << " " << base << " " << mode << endl;
        exit(EXIT_FAILURE);
    }

    return string();
} // addBase()

/**
 * Returns a suffix (mode=0) or prefix (mode=1) of @p contig with length @p length.
 *
 *	@return The suffix or prefix.
 */
string SVLocusAssembler::getEnd(const string& contig,
                                  const unsigned length,
                                  const unsigned mode) {

    const unsigned csize(contig.size());
    assert(length <= csize);

    switch (mode) {
    case 0:
        return contig.substr((csize-length),length);
    case 1:
        return contig.substr(0,length);
    default:
        cerr << "ERROR:: getEnd() : " << contig << " " << length << " " << mode << endl;
        exit(EXIT_FAILURE);
    }

    return string();
} // getEnd()

/**
 * Extends the seed contig (aka most frequent k-mer)
 *
 *	@return The extended contig.
 */
void SVLocusAssembler::walk(const string& seed,
                              const unsigned wordLength,
                              const str_uint_map_t& wordHash,
                              unsigned& stepsBackward,
                              string& contig) {

    // we start with the seed
    contig = seed;

    set<string> seenBefore;	// records k-mers already encountered during extension
    seenBefore.insert(contig);

    const str_uint_map_t::const_iterator whe(wordHash.end());

    // 0 => walk to the right, 1 => walk to the left
    for (unsigned mode(0); mode<2; ++mode) {
        while (true) {
            const string tmp(getEnd(contig,wordLength-1,mode));

#ifdef DEBUG_ASBL
            dbg_os << "# current contig : " << contig << " size : " << contig.size() << endl
                   << " getEnd : " << tmp << endl;
#endif

            if (seenBefore.find(tmp) != seenBefore.end()) {
#ifdef DEBUG_ASBL
                dbg_os << "Seen word " << tmp << " before on this walk, terminating" << endl;
#endif
                break;
            }

            seenBefore.insert(tmp);

            unsigned maxOcc(0);
            unsigned totalOcc(0);
            char maxBase('N');

            static const unsigned N_BASES(4);
            static const char BASES[] = "ACGT";

            for (unsigned i(0); i<N_BASES; ++i) {
                const char b(BASES[i]);
                const string newKey(addBase(tmp,b,mode));
#ifdef DEBUG_ASBL
                dbg_os << "Extending end : base " << b << " " << newKey << endl;
#endif
                const str_uint_map_t::const_iterator whi(wordHash.find(newKey));
                const unsigned val((whi == whe) ? 0 : whi->second);

                totalOcc += val;
                if (val > maxOcc) {
                    maxOcc  = val;
                    maxBase = b;
                }
#ifdef DEBUG_ASBL
                dbg_os << b << " " << newKey << " " << val << " " << endl;
#endif
            }
#ifdef DEBUG_ASBL
            dbg_os << "Winner is : " << maxBase << " with " << maxOcc << " occurrences." << endl;
#endif

            if ((maxOcc < minCoverage_) or
                (maxOcc < (1.-maxError_)* totalOcc)) {
#ifdef DEBUG_ASBL
                dbg_os << "Coverage or error rate below threshold.\n"
                       << "maxocc : " << maxOcc << " min coverage: " << asmOptions.minCoverage << "\n"
                       << "maxocc : " << maxOcc << " error threshold: " << ((1.0-asmOptions.maxError)* totalOcc) << endl;
#endif
                break;
            }

#ifdef DEBUG_ASBL
            dbg_os << "Adding base " << contig << " " << maxBase << " " << mode << endl;
#endif
            contig = addBase(contig,maxBase,mode);
#ifdef DEBUG_ASBL
            dbg_os << "New contig: " << contig << endl;
#endif
            if (mode == 1) {
                ++stepsBackward;
            }
        } // while(true)
#ifdef DEBUG_ASBL
        dbg_os << "mode change. Current mode " << mode << endl;
#endif
    } // for(mode)

    return;
} // walk()

void
SVLocusAssembler::getBreakendReads(const SVBreakend& bp,
								   ReadSeqVec& reads)
{
	/// define a new interval -/+ 50 bases around the center pos
	/// of the breakpoint
	static const pos_t regionSize(50);
	const pos_t centerPos(bp.interval.range.center_pos());
	const known_pos_range2 searchRange(std::max((centerPos-regionSize),0), (centerPos+regionSize));

	const unsigned bamCount(_bamStreams.size());
	for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
	{

		bam_streamer& bamStream(*_bamStreams[bamIndex]);

		// set bam stream to new search interval:
		bamStream.set_new_region(bp.interval.tid, searchRange.begin_pos(), searchRange.end_pos());

		while (bamStream.next())
		{
			const bam_record& bamRead(*(bamStream.get_record_ptr()));

			// add some criteria to filter for "interesting" reads here

			if ((bamRead.pos()-1) >= searchRange.end_pos()) break;
			reads.push_back(bamRead.get_bam_read());
		}
	}
}


void
SVLocusAssembler::assembleSVBreakend(const SVBreakend& bp,
                                     Assembly& as) {

	ReadSeqVec reads;
	getBreakendReads(bp,reads);

    /*unsigned iterations(0); // number of assembly iterations
    unsigned unused_reads_now(0);
    //cout << "Unused reads " << unused_reads_now << endl;
    unsigned unused_reads_prev(0);

    while (unused_reads_now>0)
    {
        ++iterations;
        for (unsigned int wl=wordLength_; wl<=maxWordLength_; wl+=2)
        {
            bool ret = buildContigs(data,wl,as,unused_reads_now);
            if (ret) break; // stop iteration at first successful assembly
        }
        //cout << "unused reads now " << unused_reads_now << " unused reads previous iteration " << unused_reads_prev << endl;
        if (unused_reads_now == unused_reads_prev)
        {
            // stop if no change in number of unused reads
            break;
        }
        unused_reads_prev=unused_reads_now;
        // stop assembling if we surpassed the max number of iterations
        if (iterations>maxAssemblyIterations_)
        {
            break;
        }
    }*/
}

bool SVLocusAssembler::buildContigs(const ReadSeqVec& reads,
                                    const unsigned wordLength,
                                    Assembly& as,
                                    unsigned& unused_reads) {

    //cout << "In SVLocusAssembler::buildContig. word length=" << wordLength << endl;
    const unsigned shadow_size(reads.size());

#ifdef DEBUG_ASBL
    //dbg_os << "----------------------------------------" << endl;
    for (unsigned i(0); i<shadow_size; ++i) {
        //const unsigned readNum(shadows[i].pos);
        const string& read(shadows[i].seq);
        //const bool used(shadows[i].used);
        //dbg_os << ">read" << readNum << "_used" << used << endl;
        dbg_os << read << endl;
    }
#endif

    // a set of read hashes; each read hash stores the starting positions of all kmers in the read
    vector<pair<unsigned,str_uint_map_t> > readHashes;
    // counts the number of occurrences for each kmer in all shadow reads
    str_uint_map_t wordHash;
    // most frequent kmer and its number of occurrences
    unsigned maxOcc = 0;
    string maxWord;

    for (unsigned i(0); i<shadow_size; ++i) {
        // skip shadow reads used in an previous iteration
        if (shadows[i].used) continue;
        // stores the index of a kmer in a read sequence
        //const unsigned readNum(shadows[i].pos);
        const string& read(shadows[i].seq);

        const unsigned readLen(read.size());
        assert(readLen>=wordLength);

        readHashes.resize(readHashes.size()+1);
        //readHashes.back().first=readNum;
        str_uint_map_t& readHash = readHashes.back().second;

        for (unsigned j(0); j<=(readLen-wordLength); ++j) {
            const string word(read.substr(j,wordLength));
            if (readHash.find(word) != readHash.end()) {
                // try again with different k-mer size
#ifdef DEBUG_ASBL
                dbg_os << "word " << word << " repeated in read " << read << "\n"
                       << "... try again with word length " << wordLength+2
                       << " and " << shadow_size << " shadow reads." << endl;
#endif
                return false;
            }

            // record (0-indexed) start point for word in read
            //cout << "Recording " << word << " at " << j << endl;
            readHash[word]=j;

            // count occurences
            ++wordHash[word];
            if (wordHash[word]>maxOcc) {
                //cout << "Setting max word to " << maxWord << " " << maxOcc << endl;
                maxOcc  = wordHash[word];
                maxWord = word;
            } // if (maxOcc)
        }
    }

    if (maxOcc < minCoverage_) {
#ifdef DEBUG_ASBL
        dbg_os << "Coverage too low : " << maxOcc << " " << minCoverage_ << endl;
#endif
        return false;
    }

#ifdef DEBUG_ASBL
    dbg_os << "Seeding kmer : " << maxWord << endl;
#endif

    // counts the number of steps taken backwards from the start kmer during the assembly.
    unsigned stepsBackward(0);

    // start initial assembly with most frequent kmer as seed
    AssembledContig ctg;
    walk(maxWord,wordLength,wordHash,stepsBackward,ctg.seq);

    // done with this now:
    wordHash.clear();

    const unsigned rh_size(readHashes.size());
    const unsigned contig_size(ctg.seq.size());

#ifdef DEBUG_ASBL
    dbg_os << "First pass assembly resulted in "
           << ctg.seq << "\n"
           << " with length " << contig_size << " assembled from " << rh_size << " reads. " << endl;
#endif

    for (unsigned i(0); i<rh_size; ++i) {
        //const unsigned readNum(readHashes[i].first);
        const str_uint_map_t& readHash(readHashes[i].second);
#ifdef DEBUG_ASBL
        /*{
            string read("NOTFOUND");
            for (unsigned j(0); j<shadow_size; ++j) {
                if (shadows[j].pos == readNum) {
                    read=shadows[j].seq;
                    break;
                }
            }
            dbg_os << "Checking readNum: " << readNum << " read : " << read << "\n"
                   << "vs " << ctg.seq << "\n"
                   << "number of kmers : " << readHash.size() << endl;
        }*/
#endif
        const str_uint_map_t::const_iterator rhe(readHash.end());

        // record number of reads containing the seeding kmer
        if (readHash.find(maxWord) != rhe) {
            //cout << "++seedReadCount" << endl;
            ++ctg.seedReadCount;
        }

        // store all reads aligning to the contig with gaps
        // (we don't count gaps at the end or beginning of the read)
        for (unsigned j(0); j<=(contig_size-wordLength); ++j) {
            const string word(ctg.seq.substr(j,wordLength));
            //cout << "Testing word " << word << " " << readNum << endl;
            //cout << "with counts : " << wordHash[word] << endl;
            const str_uint_map_t::const_iterator rhi(readHash.find(word));
            if (rhi != rhe) {
                //const int startPos(static_cast<int>(j)-static_cast<int>(rhi->second));
                //cout << "pos in read/contig : " << wordHash[word] << " " << j << endl;
                /*if (ctg.contigReads.find(readNum) == ctg.contigReads.end() ) {
                    ctg.contigReads[readNum] = startPos;
                    shadows[readNum].used = true;
                    --unused_reads;
                }*/
            }
        } // for each kmer...
    }

    // don't need this anymore:
    readHashes.clear();
    //cout << "final seeding reading count: " << ctg.seedReadCount << endl;
    as.push_back(ctg);
    return true;
}
