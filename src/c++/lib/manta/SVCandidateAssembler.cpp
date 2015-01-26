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
/// \author Chris Saunders
///


#include "blt_util/log.hh"
#include "blt_util/seq_util.hh"
#include "htsapi/align_path_bam_util.hh"
#include "htsapi/SimpleAlignment_bam_util.hh"
#include "manta/ShadowReadFinder.hh"
#include "manta/SVCandidateAssembler.hh"
#include "manta/SVLocusScannerSemiAligned.hh"

#include <iostream>

//#define DEBUG_REMOTES
//#define DEBUG_ASBL

SVCandidateAssembler::
SVCandidateAssembler(
    const ReadScannerOptions& scanOpt,
    const AssemblerOptions& assembleOpt,
    const AlignmentFileOptions& alignFileOpt,
    const std::string& statsFilename,
    const std::string& chromDepthFilename,
    const bam_header_info& bamHeader,
    TimeTracker& remoteTime) :
    _scanOpt(scanOpt),
    _assembleOpt(assembleOpt),
    _isAlignmentTumor(alignFileOpt.isAlignmentTumor),
    _dFilter(chromDepthFilename, scanOpt.maxDepthFactor, bamHeader),
    _dFilterRemoteReads(chromDepthFilename, scanOpt.maxDepthFactorRemoteReads, bamHeader),
    _readScanner(_scanOpt, statsFilename, alignFileOpt.alignmentFilename),
    _remoteTime(remoteTime)
{
    // setup regionless bam_streams:
    // setup all data for main analysis loop:
    for (const std::string& afile : alignFileOpt.alignmentFilename)
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
    if (! bamRead.is_paired()) return false;
    if (bamRead.is_unmapped() || bamRead.is_mate_unmapped()) return false;

    if (bamRead.map_qual() < minMapq) return false;

    if ((! isSearchForLeftOpen) && (! bamRead.is_fwd_strand())) return false;
    if ((! isSearchForRightOpen) && bamRead.is_fwd_strand()) return false;

    if (bamRead.target_id() < 0) return false;
    if (bamRead.mate_target_id() < 0) return false;

    if (bamRead.target_id() != bamRead.mate_target_id()) return true;

    /// TODO: better candidate definition based on fragment size distro:
    static const int minSize(10000);
    return (std::abs(bamRead.pos()-bamRead.mate_pos()) >= minSize);
}



/// information recorded for reads where we need to grab the mate from a remote locus
///
/// typically these are chimeras with a MAPQ0 mate used ot assemble a large insertion
///
struct RemoteReadInfo
{
    RemoteReadInfo(
        const bam_record& bamRead)
        : qname(bamRead.qname()),
          readNo(bamRead.read_no()==1 ? 2 : 1),
          tid(bamRead.mate_target_id()),
          pos(bamRead.mate_pos() - 1),
          localPos(bamRead.pos() - 1),
          readSize(bamRead.read_size()),
          isLocalFwd(bamRead.is_fwd_strand()),
          isFound(false),
          isUsed(false)
    {}

    bool
    operator<(
        const RemoteReadInfo& rhs) const
    {
        if (tid < rhs.tid) return true;
        if (tid == rhs.tid)
        {
            return (pos < rhs.pos);
        }
        return false;
    }

    std::string qname;
    int readNo; // this is read number of the target
    int tid;
    int pos;
    int localPos;
    int readSize;
    bool isLocalFwd;
    bool isFound;
    bool isUsed;
};


#ifdef DEBUG_REMOTES
static
std::ostream&
operator<<(
    std::ostream& os,
    const RemoteReadInfo& rri)
{
    os << "RRI qname: " << rri.qname
       << " readNo:" << rri.readNo
       << " tid:" << rri.tid
       << " pos:" << rri.pos
       << " localPos:" << rri.localPos
       << " readSize:" << rri.readSize
       << " isLocalFwd:" << rri.isLocalFwd
       << " isFound:" << rri.isFound
       << " isUsed:" << rri.isUsed;
    return os;
}
#endif



/// retrieve remote reads from a list of target loci in the bam
static
void
recoverRemoteReads(
    const unsigned maxNumReads,
    const bool isLocusReversed,
    const std::string& bamIndexStr,
    bam_streamer& bamStream,
    std::vector<RemoteReadInfo>& bamRemotes,
    SVCandidateAssembler::ReadIndexType& readIndex,
    AssemblyReadInput& reads,
    RemoteReadCache& remoteReadsCache)
{
    // figure out what we can handle in a single region query:
    std::sort(bamRemotes.begin(),bamRemotes.end());

    typedef std::pair<GenomeInterval, std::vector<RemoteReadInfo> > BamRegionInfo_t;
    std::vector<BamRegionInfo_t> bamRegions;

#ifdef DEBUG_REMOTES
    log_os << __FUNCTION__ << ": totalRemotes: " << bamRemotes.size() << "\n";
#endif

    int last_tid=-1;
    int last_pos=-1;
    for (const RemoteReadInfo& remote : bamRemotes)
    {
        assert(remote.tid >= 0);

        if ((last_tid == remote.tid) && (last_pos+remote.readSize >= remote.pos))
        {
            assert(! bamRegions.empty());
            GenomeInterval& interval(bamRegions.back().first);
            interval.range.set_end_pos(remote.pos);

            std::vector<RemoteReadInfo>& remotes(bamRegions.back().second);
            remotes.push_back(remote);
        }
        else
        {
            std::vector<RemoteReadInfo> remotes;
            remotes.push_back(remote);
            bamRegions.push_back(std::make_pair(GenomeInterval(remote.tid,remote.pos,remote.pos),remotes));
        }

        last_tid=remote.tid;
        last_pos=remote.pos;
    }

#ifdef DEBUG_REMOTES
    log_os << __FUNCTION__ << ": totalregions: " << bamRegions.size() << "\n";
#endif

    for (BamRegionInfo_t& bregion : bamRegions)
    {
        const GenomeInterval& interval(bregion.first);
        std::vector<RemoteReadInfo>& remotes(bregion.second);

#ifdef DEBUG_REMOTES
        log_os << __FUNCTION__ << ": begion interval " << interval << "\n";
        for (const RemoteReadInfo& remote : remotes)
        {
            log_os << " remote: " << remote.tid << " " << remote.pos << "\n";
        }

        unsigned readCount(0);
#endif

        // set bam stream to new search interval:
        bamStream.set_new_region(
            interval.tid,
            interval.range.begin_pos(),
            interval.range.end_pos()+1);

        while (bamStream.next())
        {
            if (reads.size() >= maxNumReads)
            {
#ifdef DEBUG_ASBL
                log_os << __FUNCTION__ << ": WARNING: assembly read buffer full, skipping further input\n";
#endif
                break;
            }

            const bam_record& bamRead(*(bamStream.get_record_ptr()));

            // we've gone past the last case:
            if (bamRead.pos() > (remotes.back().pos+1)) break;

            for (RemoteReadInfo& remote : remotes)
            {
#ifdef DEBUG_REMOTES
                readCount++;
                if ((readCount%1000000) == 0) log_os << " counts: " << readCount << "\n";
#endif
                if (remote.isFound) continue;
                if (bamRead.read_no() != remote.readNo) continue;
                if (bamRead.qname() != remote.qname) continue;

#ifdef DEBUG_REMOTES
                log_os << __FUNCTION__ << ": found remote: " << remote.tid << " " << remote.pos << "\n";
#endif
                remote.isFound = true;

                if (bamRead.map_qual() != 0) break;

                const char flag(bamRead.is_second() ? '2' : '1');
                const std::string readKey = remote.qname + "_" + flag + "_" + bamIndexStr;

                if (readIndex.count(readKey) != 0)
                {
#ifdef DEBUG_ASBL
                    log_os << __FUNCTION__ << ": WARNING: SmallAssembler read name collision : " << readKey << "\n";
#endif
                    break;
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

                /// add to the remote read cache used during PE scoring:
                remoteReadsCache[remote.qname] = RemoteReadPayload(bamRead.read_no(), reads.back());

                remote.isUsed = true;
                break;
            }
        }
#ifdef DEBUG_REMOTES
        log_os << __FUNCTION__ << ": total reads traversed in region: " << readCount << "\n";
#endif
    }
}



void
SVCandidateAssembler::
getBreakendReads(
    const SVBreakend& bp,
    const bool isLocusReversed,
    const reference_contig_segment& refSeq,
    const bool isSearchRemoteInsertionReads,
    RemoteReadCache& remoteReadsCache,
    ReadIndexType& readIndex,
    AssemblyReadInput& reads) const
{
    // get search range:
    known_pos_range2 searchRange;
    {
        // ideally this should be dependent on the insert size dist
        // TODO: follow-up on trial value of 200 in a separate branch/build
        // TODO: there should be a core search range and an expanded range for shadow/MAPQ0 only, shadow ranges should be left/right constrained to be consistent with center
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
    float maxDepthRemoteReads(0);
    if (isMaxDepth)
    {
        maxDepth = _dFilter.maxDepth(bp.interval.tid);
        maxDepthRemoteReads = _dFilterRemoteReads.maxDepth(bp.interval.tid);
    }
    const pos_t searchBeginPos(searchRange.begin_pos());
    const pos_t searchEndPos(searchRange.end_pos());
    std::vector<unsigned> normalDepthBuffer(searchRange.size(),0);

    bool isFirstTumor(false);

    static const unsigned maxNumReads(1000);

    const unsigned bamCount(_bamStreams.size());
    std::vector<std::vector<RemoteReadInfo> > remoteReads(bamCount);

    bool isMaxDepthRemoteReadsTriggered(false);
#ifdef FWDREV_CHECK
    /// sanity check that remote and shadow reads suggest an insertion pattern before doing an expensive remote recovery:
    std::vector<int> fwdSemiReadPos;
    std::vector<int> revSemiReadPos;
#endif

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
            if (reads.size() >= maxNumReads)
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
                if ((depthOffset >= 0) && (normalDepthBuffer[depthOffset] > maxDepthRemoteReads))
                {
                    isMaxDepthRemoteReadsTriggered=true;
                }
                if ((depthOffset >= 0) && (normalDepthBuffer[depthOffset] > maxDepth))
                {
                    continue;
                }
            }


            // check whether we do a separate search for the mate read
            if (isSearchRemoteInsertionReads)
            {
                if (isMateInsertionEvidence(bamRead, _scanOpt.minMapq, isSearchForLeftOpen, isSearchForRightOpen))
                {
#ifdef DEBUG_ASBL
                    log_os << logtag << "Adding remote bamrec: " << bamRead << '\n'
                           << "\tmapq: " << bamRead.pe_map_qual() << '\n'
                           << "\tread: " << bamRead.get_bam_read() << '\n';
#endif

                    remoteReads[bamIndex].emplace_back(bamRead);

#ifdef FWDREV_CHECK
                    if (bamRead.is_fwd_strand())
                    {
                        fwdSemiReadPos.push_back(bamRead.pos()-1);
                    }
                    else
                    {
                        revSemiReadPos.push_back(bamRead.pos()-1);
                    }
#endif
                }
            }

            // filter reads with "N"
            if (bamRead.get_bam_read().get_string().find('N') != std::string::npos) continue;

            SimpleAlignment bamAlign(getAlignment(bamRead));


            /// check for any indels in read:
            bool isIndelKeeper(false);
            if (! bamRead.is_unmapped())
            {
                using namespace ALIGNPATH;
                for (const path_segment& ps : bamAlign.path)
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

#ifdef FWDREV_CHECK
            if (isShadowKeeper)
            {
                if (bamRead.is_mate_fwd_strand())
                {
                    fwdSemiReadPos.push_back(bamRead.mate_pos()-1);
                }
                else
                {
                    revSemiReadPos.push_back(bamRead.mate_pos()-1);
                }
            }
#endif

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


    /// sanity check the remote reads to see if we're going to recover them:
    bool isRecoverRemotes(!isMaxDepthRemoteReadsTriggered);

#ifdef FWDREV_CHECK
    if (isRecoverRemotes)
    {
#ifdef DEBUG_REMOTES
        for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
        {
            unsigned fwdStrandRemotes(0);
            const std::vector<RemoteReadInfo>& bamRemotes(remoteReads[bamIndex]);
            for (const RemoteReadInfo& remote : bamRemotes)
            {
                if (remote.isLocalFwd)
                {
                    fwdStrandRemotes++;
                }
            }
            log_os << __FUNCTION__ << ": remotes for bamIndex " << bamIndex << " total: " << bamRemotes.size() << " fwd: " << fwdStrandRemotes << "\n";
        }
#endif

        // get a hack median and IQ for the remotes:
        int fwdMedian(0);
        int fwdRange(0);
        if (! fwdSemiReadPos.empty())
        {
            std::sort(fwdSemiReadPos.begin(),fwdSemiReadPos.end());
            fwdMedian=(fwdSemiReadPos[fwdSemiReadPos.size()/2]);
            fwdRange=(fwdSemiReadPos[(fwdSemiReadPos.size()*3)/4]-fwdSemiReadPos[(fwdSemiReadPos.size()*1)/4]);
        }

        int revMedian(0);
        int revRange(0);
        if (! revSemiReadPos.empty())
        {
            std::sort(revSemiReadPos.begin(),revSemiReadPos.end());
            revMedian=(revSemiReadPos[revSemiReadPos.size()/2]);
            revRange=(revSemiReadPos[(revSemiReadPos.size()*3)/4]-revSemiReadPos[(revSemiReadPos.size()*1)/4]);
        }

        if ((fwdSemiReadPos.size() <= 2) || (revSemiReadPos.size() <= 2))
        {
            isRecoverRemotes=false;
        }
        else
        {
            const int diff(revMedian-fwdMedian);
            if ((diff >= 2000) || (diff < 0))
            {
                isRecoverRemotes=false;
            }
            else if ((fwdRange >= 400) || (revRange >= 400))
            {
                isRecoverRemotes=false;
            }
        }
    }
#endif

    /// recover any remote reads:
#ifdef DEBUG_REMOTES
    log_os << __FUNCTION__ << ": isRecoverRemotes: " << isRecoverRemotes << "\n";
#endif

    if (isRecoverRemotes)
    {
        const TimeScoper remoteTIme(_remoteTime);
        for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
        {
#ifdef DEBUG_REMOTES
            log_os << __FUNCTION__ << ": starting remotes for bamindex: " << bamIndex << "\n";
#endif
            const std::string bamIndexStr(boost::lexical_cast<std::string>(bamIndex));

            bam_streamer& bamStream(*_bamStreams[bamIndex]);

            std::vector<RemoteReadInfo>& bamRemotes(remoteReads[bamIndex]);
            recoverRemoteReads(
                maxNumReads, isLocusReversed, bamIndexStr, bamStream,
                bamRemotes, readIndex, reads, remoteReadsCache);
        }
    }
}



void
SVCandidateAssembler::
assembleSingleSVBreakend(
    const SVBreakend& bp,
    const reference_contig_segment& refSeq,
    const bool isSearchRemoteInsertionReads,
    RemoteReadCache& remoteReads,
    Assembly& as) const
{
    static const bool isBpReversed(false);
    ReadIndexType readIndex;
    AssemblyReadInput reads;
    getBreakendReads(bp, isBpReversed, refSeq, isSearchRemoteInsertionReads, remoteReads, readIndex, reads);
    AssemblyReadOutput readInfo;

#ifdef ITERATIVE_ASSEMBLER
    runIterativeAssembler(_assembleOpt, reads, readInfo, as);
#else
    runSmallAssembler(_assembleOpt, reads, readInfo, as);
#endif

}



void
SVCandidateAssembler::
assembleSVBreakends(
    const SVBreakend& bp1,
    const SVBreakend& bp2,
    const bool isBp1Reversed,
    const bool isBp2Reversed,
    const reference_contig_segment& refSeq1,
    const reference_contig_segment& refSeq2,
    Assembly& as) const
{
    static const bool isSearchRemoteInsertionReads(false);
    RemoteReadCache remoteReads;
    ReadIndexType readIndex;
    AssemblyReadInput reads;
    AssemblyReadReversal readRev;
    getBreakendReads(bp1, isBp1Reversed, refSeq1, isSearchRemoteInsertionReads, remoteReads, readIndex, reads);
    readRev.resize(reads.size(),isBp1Reversed);
    getBreakendReads(bp2, isBp2Reversed, refSeq2, isSearchRemoteInsertionReads, remoteReads, readIndex, reads);
    readRev.resize(reads.size(),isBp2Reversed);
    AssemblyReadOutput readInfo;

#ifdef ITERATIVE_ASSEMBLER
    runIterativeAssembler(_assembleOpt, reads, readInfo, as);
#else
    runSmallAssembler(_assembleOpt, reads, readInfo, as);
#endif
}
