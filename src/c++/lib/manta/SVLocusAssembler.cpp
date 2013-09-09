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


#include "blt_util/align_path_bam_util.hh"
#include "blt_util/log.hh"
#include "blt_util/seq_util.hh"
#include "manta/SVLocusAssembler.hh"

#include "boost/foreach.hpp"


// compile with this macro to get verbose output:
//#define DEBUG_ASBL

#ifdef DEBUG_ASBL
#include <iostream>
std::ostream& dbg_os(std::cerr);
#endif



SVLocusAssembler::
SVLocusAssembler(
    const ReadScannerOptions& scanOpt,
    const SmallAssemblerOptions& assembleOpt,
    const AlignmentFileOptions& alignFileOpt,
    const std::string& statsFilename) :
    _scanOpt(scanOpt),
    _assembleOpt(assembleOpt),
    _readScanner(_scanOpt, statsFilename, alignFileOpt.alignmentFilename)
{
    // setup regionless bam_streams:
    // setup all data for main analysis loop:
    BOOST_FOREACH(const std::string& afile, alignFileOpt.alignmentFilename)
    {
        // avoid creating shared_ptr temporaries:
        streamPtr tmp(new bam_streamer(afile.c_str()));
        _bamStreams.push_back(tmp);
    }
}



void
SVLocusAssembler::
getBreakendReads(
    const SVBreakend& bp,
    const bool isReversed,
    ReadIndexType& readIndex,
    AssemblyReadInput& reads) const
{

    // get search range:
    known_pos_range2 searchRange;
    {
        static const size_t minIntervalSize(300);
        if (bp.interval.range.size() >= minIntervalSize)
        {
            searchRange = bp.interval.range;
        }
        else
        {
            const size_t missing = minIntervalSize - bp.interval.range.size();
            assert(missing > 0);
            const size_t wobble = missing/2;
            // FIXME : not sure what happens if (end_pos + wobble) > chromosome size?
            static const size_t zero(0);
            searchRange.set_range(std::max((bp.interval.range.begin_pos()-wobble),zero),(bp.interval.range.end_pos()+wobble));
        }
    }


    static const unsigned minClipLen(3);

    // depending on breakend type we may only be looking for candidates in one direction:
    bool isSearchForRightOpen(true);
    bool isSearchForLeftOpen(true);
    if (SVBreakendState::RIGHT_OPEN == bp.state)
    {
        isSearchForLeftOpen = false;
    }

    if (SVBreakendState::LEFT_OPEN == bp.state)
    {
        isSearchForRightOpen = false;
    }


    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const std::string bamIndexStr(boost::lexical_cast<std::string>(bamIndex));

        bam_streamer& bamStream(*_bamStreams[bamIndex]);

        // set bam stream to new search interval:
        bamStream.set_new_region(bp.interval.tid, searchRange.begin_pos(), searchRange.end_pos());

        static const unsigned MAX_NUM_READS(1000);

        //static const double minSemiAlignedScore(10.0);

        while (bamStream.next() && (reads.size() < MAX_NUM_READS))
        {
            const bam_record& bamRead(*(bamStream.get_record_ptr()));

#if 0
            // some conservative filtration criteria (including MAPQ) -- note this
            // means that if we want shadow reads this will have to be changed b/c
            // unmapped reads are filtered out here:
            if (_readScanner.isReadFiltered(bamRead)) continue;
#endif

            // include MAPQ0 because the split reads tend to have reduced mapping scores:
            if (bamRead.is_filter()) continue;
            if (bamRead.is_dup()) continue;
            if (bamRead.is_secondary()) continue;

            /// TODO: Add SA read support -- temporarily reject all supplemental reads:
            if (bamRead.is_supplement()) return;

            // FIXME: add some criteria to filter for "interesting" reads here, for now we add
            // only clipped reads and reads without N
            if ((bamRead.pos()-1) >= searchRange.end_pos()) break;

            ALIGNPATH::path_t apath;
            bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), apath);
        
            bool isClipKeeper(false);

            if (isSearchForRightOpen)
            {
                const unsigned trailingClipLen(apath_soft_clip_trail_size(apath));
                if (trailingClipLen >= minClipLen) isClipKeeper = true;
            }

            if (isSearchForLeftOpen)
            {
                const unsigned leadingClipLen(apath_soft_clip_lead_size(apath));
                if (leadingClipLen >= minClipLen) isClipKeeper = true;
            }

            //if ( !(isClipKeeper || _readScanner.isSemiAligned(bamRead) > minSemiAlignedScore )) continue;
            if (!isClipKeeper) continue;

            if (bamRead.get_bam_read().get_string().find('N') != std::string::npos) continue;

            //if ( bamRead.pe_map_qual() == 0 ) continue;
            const char flag(bamRead.is_second() ? '2' : '1');
            const std::string readKey = std::string(bamRead.qname()) + "_" + flag + "_" + bamIndexStr;

#ifdef DEBUG_ASBL
            ALIGNPATH::path_t apath;
            bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), apath);
            dbg_os << "Adding " << readKey << " " << apath << " " << bamRead.pe_map_qual() << " " << bamRead.pos() << endl;
            dbg_os << bamRead.get_bam_read().get_string() << endl;
#endif

            if (readIndex.count(readKey) == 0)
            {
                // the API gives us always the sequence w.r.t to the fwd ref, so no need
                // to reverse complement here
                readIndex.insert(std::make_pair(readKey,reads.size()));
                reads.push_back(bamRead.get_bam_read().get_string());
                if (isReversed) reverseCompStr(reads.back());
            }
            else
            {
                log_os << "WARNING : Read name collision : " << readKey << "\n";
            }
        }
    }
}



void
SVLocusAssembler::
assembleSingleSVBreakend(const SVBreakend& bp,
                         Assembly& as) const
{
    ReadIndexType readIndex;
    AssemblyReadInput reads;
    getBreakendReads(bp, false, readIndex, reads);
    AssemblyReadOutput readInfo;
    runSmallAssembler(_assembleOpt, reads, readInfo, as);
}



void
SVLocusAssembler::
assembleSVBreakends(const SVBreakend& bp1,
                    const SVBreakend& bp2,
                    const bool isBp1Reversed,
                    const bool isBp2Reversed,
                    Assembly& as) const
{
    ReadIndexType readIndex;
    AssemblyReadInput reads;
    AssemblyReadReversal readRev;
    getBreakendReads(bp1, isBp1Reversed, readIndex, reads);
    readRev.resize(reads.size(),isBp1Reversed);
    getBreakendReads(bp2, isBp2Reversed, readIndex, reads);
    readRev.resize(reads.size(),isBp2Reversed);
    AssemblyReadOutput readInfo;
    runSmallAssembler(_assembleOpt, reads, readInfo, as);
}

