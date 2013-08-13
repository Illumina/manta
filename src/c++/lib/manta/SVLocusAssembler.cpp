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

#include "blt_util/align_path_bam_util.hh"

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


SVLocusAssembler::
SVLocusAssembler(const GSCOptions& opt) :
    _scanOpt(opt.scanOpt),
    _readScanner(_scanOpt,opt.statsFilename,opt.alignmentFilename),
    // reasonable default values for 30x and 100bp reads
    _wordLength(37), _maxWordLength(61), _minContigLength(15),
    _minCoverage(1), _maxError(0.35), _minSeedReads(2),
    _maxAssemblyIterations(50)
{
    // setup regionless bam_streams:
    // setup all data for main analysis loop:
    BOOST_FOREACH(const std::string& afile, opt.alignmentFilename)
    {
        // avoid creating shared_ptr temporaries:
        streamPtr tmp(new bam_streamer(afile.c_str()));
        _bamStreams.push_back(tmp);
    }
}

/**
 * Adds base @p base to the end (mode=0) or start (mode=1) of the contig.
 *
 *	@return The extended contig.
 */
string
SVLocusAssembler::
addBase(const string& contig,
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
        cerr << "ERROR: addBase() : " << contig << " " << base << " " << mode << endl;
        exit(EXIT_FAILURE);
    }

    // FIXME: Why return empty string?
    return string();
} // addBase()

/**
 * Returns a suffix (mode=0) or prefix (mode=1) of @p contig with length @p length.
 *
 *	@return The suffix or prefix.
 */
string
SVLocusAssembler::
getEnd(const string& contig,
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
        cerr << "ERROR:: getEnd() : " << contig << " " << length << " " << mode << endl;
        exit(EXIT_FAILURE);
    }

    // FIXME: Why return empty string?
    return string();
} // getEnd()

/**
 * Extends the seed contig (aka most frequent k-mer)
 *
 */
void
SVLocusAssembler::
walk(const string& seed,
     const unsigned wordLength,
     const str_uint_map_t& wordHash,
     string& contig)
{

    // we start with the seed
    contig = seed;

    set<string> seenBefore;	// records k-mers already encountered during extension
    seenBefore.insert(contig);

    const str_uint_map_t::const_iterator whe(wordHash.end());

    // 0 => walk to the right, 1 => walk to the left
    for (unsigned mode(0); mode<2; ++mode)
    {
        while (true)
        {
            const string tmp(getEnd(contig,wordLength-1,mode));

#ifdef DEBUG_ASBL
            dbg_os << "# current contig : " << contig << " size : " << contig.size() << endl
                   << " getEnd : " << tmp << endl;
#endif

            if (seenBefore.find(tmp) != seenBefore.end())
            {
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

            for (unsigned i(0); i<N_BASES; ++i)
            {
                const char b(BASES[i]);
                const string newKey(addBase(tmp,b,mode));
#ifdef DEBUG_ASBL
                dbg_os << "Extending end : base " << b << " " << newKey << endl;
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
                dbg_os << b << " " << newKey << " " << val << " " << endl;
#endif
            }
#ifdef DEBUG_ASBL
            dbg_os << "Winner is : " << maxBase << " with " << maxOcc << " occurrences." << endl;
#endif

            if ((maxOcc < _minCoverage) or
                (maxOcc < (1.-_maxError)* totalOcc))
            {
#ifdef DEBUG_ASBL
                dbg_os << "Coverage or error rate below threshold.\n"
                       << "maxocc : " << maxOcc << " min coverage: " << _minCoverage << "\n"
                       << "maxocc : " << maxOcc << " error threshold: " << ((1.0-_maxError)* totalOcc) << endl;
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
        } // while(true)
#ifdef DEBUG_ASBL
        dbg_os << "mode change. Current mode " << mode << endl;
#endif
    } // for(mode)

    return;
} // walk()

void
SVLocusAssembler::
getBreakendReads(const SVBreakend& bp,
                 AssemblyReadMap& reads)
{
    const size_t minIntervalSize(300);
    known_pos_range2 searchRange;
    if (bp.interval.range.size() >= minIntervalSize) {
        searchRange = bp.interval.range;
    } else {
        size_t missing = minIntervalSize - bp.interval.range.size();
        assert(missing > 0);
        size_t wobble = missing/2;
        // FIXME : not sure what happens if (end_pos + wobble) > chromosome size?
        size_t zero(0);
        searchRange.set_range(std::max((bp.interval.range.begin_pos()-wobble),zero),(bp.interval.range.end_pos()+wobble));
    }

    const unsigned minClipLen(3);

    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        bam_streamer& bamStream(*_bamStreams[bamIndex]);

        // set bam stream to new search interval:
        bamStream.set_new_region(bp.interval.tid, searchRange.begin_pos(), searchRange.end_pos());

        const unsigned MAX_NUM_READS(1000);

		while (bamStream.next() && reads.size() < MAX_NUM_READS)
		{
			const bam_record& bamRead(*(bamStream.get_record_ptr()));

            // FIXME: add some criteria to filter for "interesting" reads here, for now we add
            // only clipped reads and reads without N
            if ((bamRead.pos()-1) >= searchRange.end_pos()) break;
            if (!_readScanner.isClipped(bamRead) ) continue;
            if (_readScanner.getClipLength(bamRead)<minClipLen ) continue;
            if (bamRead.get_bam_read().get_string().find('N') != std::string::npos) continue;
            //if ( bamRead.pe_map_qual() == 0 ) continue;
            string flag = "1";
            if (bamRead.is_second())
            {
                flag = "2";
            }
            string readKey = string(bamRead.qname()) + "_" + flag + "_" + boost::lexical_cast<std::string>(bamIndex);

#ifdef DEBUG_ASBL
            ALIGNPATH::path_t apath;
            bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), apath);
            dbg_os << "Adding " << readKey << " " << apath << " " << bamRead.pe_map_qual() << " " << bamRead.pos() << endl;
            dbg_os << bamRead.get_bam_read().get_string() << endl;
#endif

            if (reads.find(readKey) == reads.end()) {
                // the API gives us always the sequence w.r.t to the fwd ref, so no need
                // to reverse complement here
                reads[readKey] = AssemblyRead(bamRead.get_bam_read().get_string(),false);
            } else {
            	std::cout << "Read name collision : " << readKey << std::endl;
            }
        }
    }
}


void
SVLocusAssembler::
assembleSingleSVBreakend(const SVBreakend& bp,
                         Assembly& as)
{
    AssemblyReadMap reads;
    getBreakendReads(bp,reads);
    iterateAssembly(reads,as);
}

void
SVLocusAssembler::
assembleSVBreakends(const SVBreakend& bp1,
                    const SVBreakend& bp2,
                    Assembly& as)
{
    AssemblyReadMap reads;
    getBreakendReads(bp1,reads);
    getBreakendReads(bp2,reads);
    iterateAssembly(reads,as);
}

void
SVLocusAssembler::
iterateAssembly(AssemblyReadMap& reads, Assembly& as)
{
#ifdef DEBUG_ASBL
    std::cerr << "Starting assembly with " << reads.size() << " reads." << std::endl;
#endif

    unsigned iterations(0); // number of assembly iterations
    unsigned unused_reads_now(reads.size());
    //cout << "Unused reads " << unused_reads_now << endl;
    unsigned unused_reads_prev(reads.size());

    while (unused_reads_now>0)
    {
        ++iterations;
        for (unsigned int wl=_wordLength; wl<=_maxWordLength; wl+=2)
        {
            bool ret = buildContigs(reads,wl,as,unused_reads_now);
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
        if (iterations>_maxAssemblyIterations)
        {
#ifdef DEBUG_ASBL
            dbg_os << "Reached max number of assembly iterations: " << _maxAssemblyIterations << endl;
#endif
            break;
        }
    }
}

bool
SVLocusAssembler::
buildContigs(AssemblyReadMap& reads,
             const unsigned wordLength,
             Assembly& as,
             unsigned& unused_reads)
{

#ifdef DEBUG_ASBL
    cout << "In SVLocusAssembler::buildContig. word length=" << wordLength << endl;
    for (AssemblyReadMap::const_iterator ct = reads.begin(); ct!=reads.end(); ++ct)
    {
        dbg_os << ct->second.seq << " used=" << ct->second.used << endl;
    }
#endif

    // a set of read hashes; each read hash stores the starting positions of all kmers in the read
    vector<pair<std::string,str_uint_map_t> > readHashes;
    // counts the number of occurrences for each kmer in all shadow reads
    str_uint_map_t wordHash;
    // most frequent kmer and its number of occurrences
    unsigned maxOcc = 0;
    string maxWord;

    for (AssemblyReadMap::const_iterator ct = reads.begin(); ct != reads.end(); ++ct)
    {
        // skip shadow reads used in an previous iteration
        if (ct->second.used) continue;
        // stores the index of a kmer in a read sequence
        const string& seq(ct->second.seq);
        const string& qn(ct->first);

        const unsigned readLen(seq.size());
        assert(readLen>=wordLength);

        readHashes.resize(readHashes.size()+1);
        readHashes.back().first  = qn;
        str_uint_map_t& readHash = readHashes.back().second;

        for (unsigned j(0); j<=(readLen-wordLength); ++j)
        {
            const string word(seq.substr(j,wordLength));
            if (readHash.find(word) != readHash.end())
            {
                // try again with different k-mer size
#ifdef DEBUG_ASBL
                dbg_os << "word " << word << " repeated in read " << qn << endl
                       << "... try again with word length " << wordLength+2
                       << " and " << reads.size() << " reads left." << endl;
#endif
                return false;
            }

            // record (0-indexed) start point for word in read
            //cout << "Recording " << word << " at " << j << endl;
            readHash[word]=j;

            // count occurrences
            ++wordHash[word];
            if (wordHash[word]>maxOcc)
            {
                //cout << "Setting max word to " << maxWord << " " << maxOcc << endl;
                maxOcc  = wordHash[word];
                maxWord = word;
            } // if (maxOcc)
        }
    }

    if (maxOcc < _minCoverage)
    {
#ifdef DEBUG_ASBL
        dbg_os << "Coverage too low : " << maxOcc << " " << _minCoverage << endl;
#endif
        return false;
    }

#ifdef DEBUG_ASBL
    dbg_os << "Seeding kmer : " << maxWord << endl;
#endif

    // start initial assembly with most frequent kmer as seed
    AssembledContig ctg;
    walk(maxWord,wordLength,wordHash,ctg.seq);

    // done with this now:
    wordHash.clear();

    const unsigned rh_size(readHashes.size());
    const unsigned contig_size(ctg.seq.size());
    const AssemblyReadMap::const_iterator arm_e = reads.end();

#ifdef DEBUG_ASBL
    dbg_os << "First pass assembly resulted in "
           << ctg.seq << "\n"
           << " with length " << contig_size << ". Input consisted of " << rh_size << " reads. " << endl;
#endif

    for (unsigned i(0); i<rh_size; ++i)
    {
        const std::string qName(readHashes[i].first);
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
            const string word(ctg.seq.substr(j,wordLength));
            //cout << "Testing word " << word << " " << readNum << endl;
            //cout << "with counts : " << wordHash[word] << endl;
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
                    std::cerr << "Read " << qName << " not found in hash" << std::endl;
                }
            }
        } // for each kmer...
    }

    // don't need this anymore:
    readHashes.clear();
#ifdef DEBUG_ASBL
    cout << "final seeding reading count: " << ctg.seedReadCount << endl;
#endif
    if (ctg.seedReadCount > _minSeedReads)
    {
        as.push_back(ctg);
    }
    else
    {
#ifdef DEBUG_ASBL
        cout << "which is below minSeedReadCount of " << _minSeedReads << " discarding."  << endl;
#endif
    }
    return true;
}
