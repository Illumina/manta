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


#include "blt_util/align_path_bam_util.hh"
#include "blt_util/log.hh"
#include "blt_util/seq_util.hh"
#include "manta/ShadowReadFinder.hh"
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



static
bool
isMateInsertionEvidence(
    const bam_record& bamRead,
    const unsigned minMapq,
    const bool isSearchForLeftOpen,
    const bool isSearchForRightOpen)
{
    if (bamRead.is_unmapped() || bamRead.is_mate_unmapped()) return false;

    if (bamRead.map_qual() < minMapq) return false;

    if ((! isSearchForLeftOpen) && (! bamRead.is_fwd_strand())) return false;
    if ((! isSearchForRightOpen) && bamRead.is_fwd_strand()) return false;

    if (bamRead.target_id() != bamRead.mate_target_id()) return true;

    /// TODO: better candidate definition based on fragment size distro:
    static const int minSize(10000);
    return (std::abs(bamRead.pos()-bamRead.mate_pos()) >= minSize);
}



struct RemoteReadInfo
{
    RemoteReadInfo(
        const bam_record& bamRead)
      : qname(bamRead.qname()),
        readNo(bamRead.read_no()==1 ? 2 : 1),
        tid(bamRead.mate_target_id()),
        pos(bamRead.mate_pos() - 1)
    {}

    std::string qname;
    int readNo; // this is read number of the target
    int tid;
    int pos;
};



void
SVLocusAssembler::
getBreakendReads(
    const SVBreakend& bp,
    const bool isLocusReversed,
    const reference_contig_segment& refSeq,
    const bool isSearchRemoteInsertionReads,
    ReadIndexType& readIndex,
    AssemblyReadInput& reads) const
{
    // get search range:
    known_pos_range2 searchRange;
    {
        // ideally this should be dependent on the insert size dist
        // TODO: follow-up on trial value of 200 in a separate branch/build
        // TODO: there should be a core search range and and expanded range for shadow/MAPQ0 only, shadow ranges should be left/right constrained to be consistent with center
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
    float maxDepth(0);
    if (isMaxDepth)
    {
        maxDepth = _dFilter.maxDepth(bp.interval.tid);
    }
    const pos_t searchBeginPos(searchRange.begin_pos());
    const pos_t searchEndPos(searchRange.end_pos());
    std::vector<unsigned> normalDepthBuffer(searchRange.size(),0);

    bool isFirstTumor(false);

    static const unsigned MAX_NUM_READS(1000);

    const unsigned bamCount(_bamStreams.size());
    std::vector<std::vector<RemoteReadInfo> > remoteReads(bamCount);

    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const bool isTumor(_isAlignmentTumor[bamIndex]);

        // assert that the expected sample order is all normal samples first,
        // followed by all tumor samples
        if (isTumor) isFirstTumor=true;
        assert((! isFirstTumor) || isTumor);

        const std::string bamIndexStr(boost::lexical_cast<std::string>(bamIndex));

        bam_streamer& bamStream(*_bamStreams[bamIndex]);

        // set bam stream to new search interval:
        bamStream.set_new_region(bp.interval.tid, searchBeginPos, searchEndPos);

        ShadowReadFinder shadow(_scanOpt.minSingletonMapqCandidates, isSearchForLeftOpen, isSearchForRightOpen);

#ifdef DEBUG_ASBL
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

            // don't filter out MAPQ0 because the split reads tend to have reduced mapping scores:
            if (SVLocusScanner::isReadFilteredCore(bamRead)) continue;

            const pos_t refPos(bamRead.pos()-1);
            if (refPos >= searchEndPos) break;

            if (isMaxDepth)
            {
                assert(refPos<searchEndPos);
                const pos_t depthOffset(refPos - searchBeginPos);
                if ((depthOffset >= 0) && (normalDepthBuffer[depthOffset] > maxDepth)) continue;
            }


            // check whether we do a separate search for the mate read
            if (isSearchRemoteInsertionReads)
            {
                if (isMateInsertionEvidence(bamRead, _scanOpt.minMapq, isSearchForLeftOpen, isSearchForRightOpen))
                {
                    remoteReads[bamIndex].push_back(RemoteReadInfo(bamRead));
                }
            }

            // filter reads with "N"
            if (bamRead.get_bam_read().get_string().find('N') != std::string::npos) continue;

            SimpleAlignment bamAlign(bamRead);


            /// check for any indels in read:
            bool isIndelKeeper(false);
            if (! bamRead.is_unmapped())
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

            /// this test covered semi-aligned and soft-clip and split reads together
            bool isSemiAlignedKeeper(false);
            if (! bamRead.is_unmapped())
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

            const bool isShadowKeeper(shadow.check(bamRead));

            if (! (isIndelKeeper
                   || isSemiAlignedKeeper
                   || isShadowKeeper
                  )) continue;

#ifdef DEBUG_ASBL
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
            log_os << "isIndelKeeper: " << isIndelKeeper
                   << " isSemiAlignedKeeper: " << isSemiAlignedKeeper
                   << " isShadowKeeper: " << isShadowKeeper
                   << '\n';
#endif

            if (readIndex.count(readKey) != 0)
            {
#ifdef DEBUG_ASBL
                log_os << logtag << "WARNING: SmallAssembler read name collision : " << readKey << "\n";
#endif
                continue;
            }

            readIndex.insert(std::make_pair(readKey,reads.size()));
            reads.push_back(bamRead.get_bam_read().get_string());

            bool isReversed(isLocusReversed);

            // if shadow read, determine if we need to reverse:
            if (isShadowKeeper)
            {
                if (bamRead.is_mate_fwd_strand())
                {
                    isReversed = (! isReversed);
                }
            }

            if (isReversed) reverseCompStr(reads.back());
        }

#ifdef DEBUG_ASBL
        log_os << logtag << "bam " << bamIndex
               << " indel: " << indelCount
               << " semi-aligned " << semiAlignedCount
               << " shadow " << shadowCount
               << '\n';
#endif
    }


    /// recover any remote reads:
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const std::string bamIndexStr(boost::lexical_cast<std::string>(bamIndex));

        bam_streamer& bamStream(*_bamStreams[bamIndex]);

        const std::vector<RemoteReadInfo>& bamRemotes(remoteReads[bamIndex]);

        BOOST_FOREACH(const RemoteReadInfo& remote, bamRemotes)
        {
            // set bam stream to new search interval:
            bamStream.set_new_region(remote.tid, remote.pos-1, remote.pos+1);

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

                if (bamRead.qname() != remote.qname) continue;
                if (bamRead.read_no() != remote.readNo) continue;
                if (bamRead.map_qual() != 0) break;

                const char flag(bamRead.is_second() ? '2' : '1');
                const std::string readKey = remote.qname + "_" + flag + "_" + bamIndexStr;

                if (readIndex.count(readKey) != 0)
                {
    #ifdef DEBUG_ASBL
                    log_os << logtag << "WARNING: SmallAssembler read name collision : " << readKey << "\n";
    #endif
                    continue;
                }

                readIndex.insert(std::make_pair(readKey,reads.size()));
                reads.push_back(bamRead.get_bam_read().get_string());

                bool isReversed(isLocusReversed);

                // determine if we need to reverse:
                if (bamRead.is_fwd_strand() == bamRead.is_mate_fwd_strand())
                {
                    isReversed = (! isReversed);
                }

                if (isReversed) reverseCompStr(reads.back());
                break;
            }
        }
    }
}



void
SVLocusAssembler::
assembleSingleSVBreakend(
    const SVBreakend& bp,
    const reference_contig_segment& refSeq,
    const bool isSearchRemoteInsertionReads,
    Assembly& as) const
{
    static const bool isBpReversed(false);
    ReadIndexType readIndex;
    AssemblyReadInput reads;
    getBreakendReads(bp, isBpReversed, refSeq, isSearchRemoteInsertionReads, readIndex, reads);
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
    static const bool isSearchRemoteInsertionReads(false);
    ReadIndexType readIndex;
    AssemblyReadInput reads;
    AssemblyReadReversal readRev;
    getBreakendReads(bp1, isBp1Reversed, refSeq1, isSearchRemoteInsertionReads, readIndex, reads);
    readRev.resize(reads.size(),isBp1Reversed);
    getBreakendReads(bp2, isBp2Reversed, refSeq2, isSearchRemoteInsertionReads, readIndex, reads);
    readRev.resize(reads.size(),isBp2Reversed);
    AssemblyReadOutput readInfo;
    runSmallAssembler(_assembleOpt, reads, readInfo, as);
}
