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
SVScorer(
    const GSCOptions& opt,
    const bam_header_info& header) :
    _isAlignmentTumor(opt.alignFileOpt.isAlignmentTumor),
    _somaticOpt(opt.somaticOpt),
    _dFilter(opt.chromDepthFilename, opt.somaticOpt.maxDepthFactor, header),
    _readScanner(opt.scanOpt,opt.statsFilename,opt.alignFileOpt.alignmentFilename)
{
    // setup regionless bam_streams:
    // setup all data for main analysis loop:
    BOOST_FOREACH(const std::string& afile, opt.alignFileOpt.alignmentFilename)
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
scoreSplitReads(
		const SVBreakend& bp,
		const AssembledContig& contig,
		const unsigned bp1ContigOfs,
		const unsigned bp2ContigOfs,
		const reference_contig_segment& bp1Ref,
		const reference_contig_segment& bp2Ref,
		const unsigned bp1RefOfs,
		const unsigned bp2RefOfs,
		bam_streamer& readStream,
		const bool isTumor,
		SomaticSVScoreInfo& ssInfo)
{

	// extract reads overlapping the break point
	readStream.set_new_region(bp.interval.tid, bp.interval.range.begin_pos(), bp.interval.range.end_pos());
	while (readStream.next())
	{
		const bam_record& bamRead(*(readStream.get_record_ptr()));
		std::string readSeq = bamRead.get_bam_read().get_string();
		splitReadAlignment bp1ContigAlign = alignSplitRead(readSeq, contig, bp1ContigOfs);
		splitReadAlignment bp2ContigAlign = alignSplitRead(readSeq, contig, bp2ContigOfs);
		splitReadAlignment bp1RefAlign = alignSplitRead(readSeq, bp1Ref, bp1RefOfs);
		splitReadAlignment bp2RefAlign = alignSplitRead(readSeq, bp2Ref, bp2RefOfs);

		// scoring
		float bp1ContigScore(0);
		float bp2ContigScore(0);
		float bp1RefScore(0);
		float bp2RefScore(0);

		if (bp1ContigAlign.has_evidence())
			bp1ContigScore = bp1ContigAlign.get_score();
		if (bp2ContigAlign.has_evidence())
			bp2ContigScore = bp2ContigAlign.get_score();

		if (bp1RefAlign.has_evidence())
			bp1RefScore = bp1RefAlign.get_score();
		if (bp2RefAlign.has_evidence())
			bp2RefScore = bp2RefAlign.get_score();

		if (isTumor)
		{
			ssInfo.tumor.contigSplitReads += std::max(bp1ContigScore, bp2ContigScore);
			ssInfo.tumor.refSplitReads += std::max(bp1RefScore, bp2RefScore);
		}
		else
		{
			ssInfo.normal.contigSplitReads += std::max(bp1ContigScore, bp2ContigScore);
			ssInfo.normal.refSplitReads += std::max(bp1RefScore, bp2RefScore);
		}
	}
}


void
SVScorer::
scoreSomaticSV(
    const SVCandidateSetData& svData,
    const SVCandidateAssemblyData& assemblyData,
    const SVCandidate& sv,
    SomaticSVScoreInfo& ssInfo)
{
    ssInfo.clear();

    // get somatic contig
    const AssembledContig& contig = assemblyData.contigs[assemblyData.bestAlignmentIndex];
    const JumpAlignmentResult<int>& alignment = assemblyData.spanningAlignments[assemblyData.bestAlignmentIndex];
    // get offsets of breakpoints in the contig
    const unsigned align1Size(apath_read_length(alignment.align1.apath));
    const unsigned insertSize(alignment.jumpInsertSize);
    const unsigned bp1ContigOffset = align1Size - 1;
    const unsigned bp2ContigOffset = align1Size + insertSize - 1;

    // get reference regions
    const reference_contig_segment& bp1Ref = assemblyData.bp1ref;
    const reference_contig_segment& bp2Ref = assemblyData.bp2ref;
    // get offsets of breakpoints in the reference regions
    unsigned bp1RefOffset = sv.bp1.interval.range.begin_pos() - bp1Ref.get_offset();
    unsigned bp2RefOffset = sv.bp2.interval.range.begin_pos() - bp2Ref.get_offset();

    // first exercise -- just count the sample assignment of the pairs we already have:
    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const bool isTumor(_isAlignmentTumor[bamIndex]);

        // consider 2-locus events first
        // TODO: to add local assembly later
        if (assemblyData.isSpanning)
        {
			streamPtr& bamPtr(_bamStreams[bamIndex]);
			bam_streamer& read_stream(*bamPtr);
			// scoring split reads overlapping bp1
			scoreSplitReads(sv.bp1, contig, bp1ContigOffset, bp2ContigOffset,
					        bp1Ref, bp2Ref, bp1RefOffset, bp2RefOffset,
					        read_stream, isTumor, ssInfo);
			// scoring split reads overlapping bp2
			scoreSplitReads(sv.bp2, contig, bp1ContigOffset, bp2ContigOffset,
			        		bp1Ref, bp2Ref, bp1RefOffset, bp2RefOffset,
			        		read_stream, isTumor, ssInfo);
        }

        SVSampleInfo& sample(isTumor ? ssInfo.tumor : ssInfo.normal);
        const SVCandidateSetReadPairSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
        BOOST_FOREACH(const SVCandidateSetReadPair& pair, svDataGroup)
        {
            if (0 == std::count(pair.svIndex.begin(),pair.svIndex.end(), sv.candidateIndex)) continue;

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
    ssInfo.bp1MaxDepth=(getBreakendMaxMappedDepth(sv.bp1));
    ssInfo.bp2MaxDepth=(getBreakendMaxMappedDepth(sv.bp2));

    // Get Data on standard read pairs crossing the two breakends,
    // and get a breakend depth estimate

    // apply filters
    if (_dFilter.isMaxDepthFilter())
    {
        if (ssInfo.bp1MaxDepth > _dFilter.maxDepth(sv.bp1.interval.tid))
        {
            ssInfo.filters.insert(_somaticOpt.maxDepthFilterLabel);
        }
        else if (ssInfo.bp2MaxDepth > _dFilter.maxDepth(sv.bp2.interval.tid))
        {
            ssInfo.filters.insert(_somaticOpt.maxDepthFilterLabel);
        }

    }

    // assign bogus somatic score just to get started:
    bool isSomatic(true);
    if ((ssInfo.normal.spanPairs > 1))
    	isSomatic=false;

    if (isSomatic)
    {
        const bool lowPairSupport(ssInfo.tumor.spanPairs < 6);
        const bool lowSingleSupport((ssInfo.tumor.bp1SpanReads < 14) || (ssInfo.tumor.bp2SpanReads < 14));
        const bool highSingleContam((ssInfo.normal.bp1SpanReads > 1) || (ssInfo.normal.bp2SpanReads > 1));

        /// allow single pair support to rescue an SV only if the evidence looks REALLY good:
        if (lowPairSupport && (lowSingleSupport || highSingleContam))
        	isSomatic=false;
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




