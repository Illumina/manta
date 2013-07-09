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
/// \author Chris Saunders
///

#include "SVScorer.hh"

#include "blt_util/align_path_bam_util.hh"
#include "blt_util/bam_streamer.hh"
#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "manta/ReadGroupStatsSet.hh"

#include "boost/foreach.hpp"

#include <algorithm>
#include <iostream>



SVScorer::
SVScorer(const GSCOptions& opt) :
    _isAlignmentTumor(opt.isAlignmentTumor),
    _somaticOpt(opt.somaticOpt),
    _readScanner(opt.scanOpt,opt.statsFilename,opt.alignmentFilename)
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



static
void
addReadToDepthEst(
    const bam_record& bamRead,
    const pos_t beginPos,
    std::vector<unsigned>& depth)
{
    using namespace ALIGNPATH;

    const pos_t endPos(beginPos+depth.size());

    // get cigar:
    path_t apath;
    bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), apath);

    pos_t refPos(bamRead.pos()-1);
    BOOST_FOREACH(const path_segment& ps, apath)
    {
        if (refPos>=endPos) return;

        if (MATCH == ps.type)
        {
            for (pos_t pos(refPos); pos < (refPos+static_cast<pos_t>(ps.length)); ++pos)
            {
                if (pos>=beginPos)
                {
                    if (pos>=endPos) return;
                    depth[pos-beginPos]++;
                }
            }
        }
        if (is_segment_type_ref_length(ps.type)) refPos += ps.length;
    }
}



unsigned
SVScorer::
getBreakendMaxMappedDepth(const SVBreakend& bp)
{
    /// define a new interval -/+ 50 bases around the center pos
    /// of the breakpoint
    static const pos_t regionSize(50);
    const pos_t centerPos(bp.interval.range.center_pos());
    const known_pos_range2 searchRange(std::max((centerPos-regionSize),0), (centerPos+regionSize));

    std::vector<unsigned> depth(searchRange.size(),0);

    bool isNormalFound(false);

    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        if (_isAlignmentTumor[bamIndex]) continue;

        bam_streamer& bamStream(*_bamStreams[bamIndex]);

        // set bam stream to new search interval:
        bamStream.set_new_region(bp.interval.tid, searchRange.begin_pos(), searchRange.end_pos());

        while (bamStream.next())
        {
            const bam_record& bamRead(*(bamStream.get_record_ptr()));

            // turn filtration off down to mapped only to match depth estimate method:
            //if (_readScanner.isReadFiltered(bamRead)) continue;
            if (bamRead.is_unmapped()) continue;

            if ((bamRead.pos()-1) >= searchRange.end_pos()) break;

            addReadToDepthEst(bamRead,searchRange.begin_pos(),depth);
        }

        isNormalFound=true;
        break;
    }

    assert(isNormalFound);

    return *(std::max_element(depth.begin(),depth.end()));
}




void
SVScorer::
scoreSomaticSV(
    const SVCandidateData& svData,
    const unsigned svIndex,
    const SVCandidate& sv,
    SomaticSVScoreInfo& ssInfo)
{
    ssInfo.clear();

    // first exercise -- just count the sample assignment of the pairs we already have:

    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const bool isTumor(_isAlignmentTumor[bamIndex]);
        SVSampleInfo& sample(isTumor ? ssInfo.tumor : ssInfo.normal);

        const SVCandidateDataGroup& svDataGroup(svData.getDataGroup(bamIndex));
        BOOST_FOREACH(const SVCandidateReadPair& pair, svDataGroup)
        {
            if (svIndex != pair.svIndex) continue;

            if (pair.read1.isSet())
            {
                sample.bp1SpanReads += 1;
            }
            if (pair.read2.isSet())
            {
                sample.bp2SpanReads += 1;
            }
            if (pair.read1.isSet() && pair.read2.isSet())
            {
                sample.spanPairs += 1;
            }
        }
    }

    // get breakend center_pos depth estimate:
    const unsigned bp1MaxDepth(getBreakendMaxMappedDepth(sv.bp1));
    const unsigned bp2MaxDepth(getBreakendMaxMappedDepth(sv.bp2));

    if (bp1MaxDepth+bp2MaxDepth > 1) return;

    // Get Data on standard read pairs crossing the two breakends,
    // and get a breakend depth estimate

    // apply filters


    // assign bogus somatic score just to get started:
    bool isSomatic(true);
    if (ssInfo.normal.spanPairs > 1) isSomatic=false;

    if (isSomatic)
    {
        const bool lowPairSupport(ssInfo.tumor.spanPairs < 6);
        const bool lowSingleSupport((ssInfo.tumor.bp1SpanReads < 14) || (ssInfo.tumor.bp2SpanReads < 14));
        const bool highSingleContam((ssInfo.normal.bp1SpanReads > 1) || (ssInfo.normal.bp2SpanReads > 1));

        /// allow single pair support to rescue an SV only if the evidence looks REALLY good:
        if (lowPairSupport && (lowSingleSupport || highSingleContam)) isSomatic=false;
    }

    if (isSomatic)
    {
        if (ssInfo.normal.spanPairs)
        {
            const double ratio(static_cast<double>(ssInfo.tumor.spanPairs)/static_cast<double>(ssInfo.normal.spanPairs));
            if (ratio<9)
            {
                isSomatic=false;
            }
        }
        if (ssInfo.normal.bp1SpanReads)
        {
            const double ratio(static_cast<double>(ssInfo.tumor.bp1SpanReads)/static_cast<double>(ssInfo.normal.bp1SpanReads));
            if (ratio<9)
            {
                isSomatic=false;
            }
        }
        if (ssInfo.normal.bp2SpanReads)
        {
            const double ratio(static_cast<double>(ssInfo.tumor.bp2SpanReads)/static_cast<double>(ssInfo.normal.bp2SpanReads));
            if (ratio<9)
            {
                isSomatic=false;
            }
        }
    }

    if (isSomatic) ssInfo.somaticScore=60;
}




