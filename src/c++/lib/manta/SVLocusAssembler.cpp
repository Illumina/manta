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
// <https://github.com/sequencing/licenses/>
//

///
/// \author Ole Schulz-Trieglaff
///


#include "blt_util/align_path_bam_util.hh"
#include "blt_util/log.hh"
#include "blt_util/seq_util.hh"
#include "manta/SVLocusAssembler.hh"
#include "manta/SVLocusScannerSemiAligned.hh"

#include "boost/foreach.hpp"

#include <iostream>


//#define DEBUG_ASBL



SVLocusAssembler::
SVLocusAssembler(
    const ReadScannerOptions& scanOpt,
    const SmallAssemblerOptions& assembleOpt,
    const AlignmentFileOptions& alignFileOpt,
    const std::string& statsFilename,
    const std::string& chromDepthFilename,
    const bam_header_info& bamHeader) :
    _scanOpt(scanOpt),
    _assembleOpt(assembleOpt),
    _isAlignmentTumor(alignFileOpt.isAlignmentTumor),
    _dFilter(chromDepthFilename, scanOpt.maxDepthFactor, bamHeader),
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



/// approximate depth tracking -- don't bother reading the cigar string, just assume a perfect match of
/// size read_size
static
void
addReadToDepthEst(
    const bam_record& bamRead,
    const pos_t beginPos,
    std::vector<unsigned>& depth)
{
    const pos_t endPos(beginPos+depth.size());
    const pos_t refStart(bamRead.pos()-1);

    const pos_t readSize(bamRead.read_size());
    for (pos_t readIndex(std::max(0,(beginPos-refStart))); readIndex<readSize; ++readIndex)
    {
        const pos_t refPos(refStart+readIndex);
        if (refPos>=endPos) return;
        const pos_t depthIndex(refPos-beginPos);
        assert(depthIndex>=0);

        depth[depthIndex]++;
    }
}



void
SVLocusAssembler::
getBreakendReads(
    const SVBreakend& bp,
    const bool isReversed,
    const reference_contig_segment& refSeq,
    ReadIndexType& readIndex,
    AssemblyReadInput& reads) const
{
    // get search range:
    known_pos_range2 searchRange;
    {
        // ideally this should be dependent on the insert size dist
        // TODO: follow-up on trial value of 200 in a separate branch/build
        static const size_t minIntervalSize(400);
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

#ifdef DEBUG_ASBL
    static const std::string logtag("SVLocusAssembler::getBreakendReads: ");
    log_os << logtag << "searchRange " << searchRange << "\n";
#endif

    // for assembler reads, look for indels at report size or somewhat smaller
    const unsigned minAssembleIndelSize(_scanOpt.minCandidateVariantSize/2);

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

    const bool isMaxDepth(_dFilter.isMaxDepthFilter());
    const float maxDepth(_dFilter.maxDepth(bp.interval.tid));
    const pos_t searchBeginPos(searchRange.begin_pos());
    const pos_t searchEndPos(searchRange.end_pos());
    std::vector<unsigned> normalDepthBuffer(searchRange.size(),0);

    bool isFirstTumor(false);

    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const bool isTumor(_isAlignmentTumor[bamIndex]);

        /// assert expected sample order of all normal, then all tumor:
        if (isTumor) isFirstTumor=true;
        assert((! isFirstTumor) || isTumor);

        const std::string bamIndexStr(boost::lexical_cast<std::string>(bamIndex));

        bam_streamer& bamStream(*_bamStreams[bamIndex]);

        // set bam stream to new search interval:
        bamStream.set_new_region(bp.interval.tid, searchBeginPos, searchEndPos);

        // Singleton/shadow pairs *should* appear consecutively in the BAM file
        // this keeps track of the mapq score of the singleton read such that
        // we can access it when we look at the shadow
        bool isLastSet(false);
        uint8_t lastMapq(0);
        std::string lastQname;

        static const unsigned MAX_NUM_READS(1000);

#ifdef DEBUG_ASBL
        unsigned clipCount(0);
        unsigned indelCount(0);
        unsigned semiAlignedCount(0);
        unsigned shadowCount(0);
#endif

        while (bamStream.next())
        {
            if (reads.size() >= MAX_NUM_READS)
            {
#ifdef DEBUG_ASBL
                log_os << logtag << "WARNING: assembly read buffer full, skipping further input\n";
#endif
                break;
            }

            const bam_record& bamRead(*(bamStream.get_record_ptr()));

            if (isMaxDepth)
            {
                if (! isTumor)
                {
                    // depth estimation relies on a simple filtration criteria to stay in sync with the chromosome mean
                    // depth estimates:
                    if (! bamRead.is_unmapped())
                    {
                        addReadToDepthEst(bamRead, searchBeginPos, normalDepthBuffer);
                    }
                }
            }

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
            if (bamRead.is_supplement()) continue;

            const pos_t refPos(bamRead.pos()-1);
            if (refPos >= searchEndPos) break;

            if (isMaxDepth)
            {
                assert(refPos<searchEndPos);
                const pos_t depthOffset(refPos - searchBeginPos);
                if ((depthOffset >= 0) && (normalDepthBuffer[depthOffset] > maxDepth)) continue;
            }

            // filter reads with "N"
            if (bamRead.get_bam_read().get_string().find('N') != std::string::npos) continue;

            SimpleAlignment bamAlign(bamRead);

            /// check whether we keep this read because of soft clipping:
            bool isClipKeeper(false);
            {
                static const unsigned minSoftClipLen(4);

                unsigned leadingClipLen(0);
                unsigned trailingClipLen(0);
                getSVBreakendCandidateClip(bamRead, bamAlign.path, leadingClipLen, trailingClipLen);

                if (isSearchForRightOpen)
                {
                    if (trailingClipLen >= minSoftClipLen) isClipKeeper = true;
                }

                if (isSearchForLeftOpen)
                {
                    if (leadingClipLen >= minSoftClipLen) isClipKeeper = true;
                }
            }

            /// check for any indels in read:
            bool isIndelKeeper(false);
            {
                using namespace ALIGNPATH;
                BOOST_FOREACH(const path_segment& ps, bamAlign.path)
                {
                    if (is_segment_type_indel(ps.type))
                    {
                        if (ps.length>=minAssembleIndelSize) isIndelKeeper = true;
                        break;
                    }
                }
            }


            bool isSemiAlignedKeeper(false);
            {
                static const unsigned minMismatchLen(4);

                unsigned leadingMismatchLen(0);
                unsigned trailingMismatchLen(0);
                getSVBreakendCandidateSemiAligned(bamRead, bamAlign, refSeq, leadingMismatchLen, trailingMismatchLen);

                if (isSearchForRightOpen)
                {
                    if (trailingMismatchLen >= minMismatchLen) isSemiAlignedKeeper = true;
                }

                if (isSearchForLeftOpen)
                {
                    if (leadingMismatchLen >= minMismatchLen) isSemiAlignedKeeper = true;
                }

#if 0
                if (isSemiAligned(bamRead,ref,_scanOpt.minSemiAlignedScoreCandidates))
                {
                    isSemiAlignedKeeper = true;
                }
#endif
            }

            bool isShadowKeeper(false);

            if (isLastSet)
            {
                if (isGoodShadow(bamRead,
                                 lastMapq,
                                 lastQname,
                                 _scanOpt.minSingletonMapqCandidates))
                {
                    isShadowKeeper = true;
                }
            }


            lastMapq  = bamRead.map_qual();
            lastQname = bamRead.qname();
            isLastSet = true;

            if (! (isClipKeeper
                   || isIndelKeeper
                   || isSemiAlignedKeeper
                   || isShadowKeeper
                  )) continue;

#ifdef DEBUG_ASBL
            if (isClipKeeper) ++clipCount;
            if (isIndelKeeper) ++indelCount;
            if (isSemiAlignedKeeper) ++semiAlignedCount;
            if (isShadowKeeper) ++shadowCount;
#endif
            //if ( bamRead.pe_map_qual() == 0 ) continue;
            const char flag(bamRead.is_second() ? '2' : '1');
            const std::string readKey = std::string(bamRead.qname()) + "_" + flag + "_" + bamIndexStr;

#ifdef DEBUG_ASBL
            log_os << logtag << "Adding bamrec: " << bamRead << '\n'
                   << "\tmapq: " << bamRead.pe_map_qual() << '\n'
                   << "\tread: " << bamRead.get_bam_read() << '\n';
            log_os << "isClipKeeper: " << isClipKeeper
                   << " isIndelKeeper: " << isIndelKeeper
                   << " isSemiAlignedKeeper: " << isSemiAlignedKeeper
                   << " isShadowKeeper: " << isShadowKeeper
                   << '\n';
#endif

            if (readIndex.count(readKey) == 0)
            {
                readIndex.insert(std::make_pair(readKey,reads.size()));
                reads.push_back(bamRead.get_bam_read().get_string());
                if (isReversed) reverseCompStr(reads.back());
            }
            else
            {
#ifdef DEBUG_ASBL
                log_os << logtag << "WARNING: SmallAssembler read name collision : " << readKey << "\n";
#endif
            }
        }

#ifdef DEBUG_ASBL
        log_os << logtag << "bam " << bamIndex
               << " clip: " << clipCount
               << " indel: " << indelCount
               << " semi-aligned " << semiAlignedCount
               << " shadow " << shadowCount
               << '\n';
#endif
    }
}



void
SVLocusAssembler::
assembleSingleSVBreakend(
    const SVBreakend& bp,
    const reference_contig_segment& refSeq,
    Assembly& as) const
{
    static const bool isBpReversed(false);
    ReadIndexType readIndex;
    AssemblyReadInput reads;
    getBreakendReads(bp, isBpReversed, refSeq, readIndex, reads);
    AssemblyReadOutput readInfo;
    runSmallAssembler(_assembleOpt, reads, readInfo, as);
}



void
SVLocusAssembler::
assembleSVBreakends(const SVBreakend& bp1,
                    const SVBreakend& bp2,
                    const bool isBp1Reversed,
                    const bool isBp2Reversed,
                    const reference_contig_segment& refSeq1,
                    const reference_contig_segment& refSeq2,
                    Assembly& as) const
{
    ReadIndexType readIndex;
    AssemblyReadInput reads;
    AssemblyReadReversal readRev;
    getBreakendReads(bp1, isBp1Reversed, refSeq1, readIndex, reads);
    readRev.resize(reads.size(),isBp1Reversed);
    getBreakendReads(bp2, isBp2Reversed, refSeq2, readIndex, reads);
    readRev.resize(reads.size(),isBp2Reversed);
    AssemblyReadOutput readInfo;
    runSmallAssembler(_assembleOpt, reads, readInfo, as);
}

