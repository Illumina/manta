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
/// \author Chris Saunders and Xiaoyu Chen
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
#include <string>

#define DEBUG_SVS


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
		bool isBp1,
		const SVBreakend& bp,
		const SVAlignmentInfo& svAlignInfo,
		read_map_t& readMap,
		bam_streamer& readStream,
		SVSampleInfo& sample)
{
	// extract reads overlapping the break point
	readStream.set_new_region(bp.interval.tid, bp.interval.range.begin_pos(), bp.interval.range.end_pos());
	while (readStream.next())
	{
		const bam_record& bamRead(*(readStream.get_record_ptr()));
		const std::string readSeq = bamRead.get_bam_read().get_string();
		const std::string readSeqRC = bamRead.get_bam_read().get_rc_string();
		const unsigned readMapQ = bamRead.map_qual();

		// TODO: the logic of R1 & R2 should change for alignments with split-reads (e.g. BWA-MEME)
		std::string readId(bamRead.qname());
		if (bamRead.is_first())
			readId += "_R1";
		else if(bamRead.is_second())
			readId += "_R2";

		// map all reads for bp1
		// then for bp2, skip reads that have been considered for bp1 to avoid double-counting
		if (isBp1)
		{
			readMap[readId] = true;
			log_os << readId << "mapped.\n";
		}
		else if (readMap.find(readId) != readMap.end())
		{
			log_os << readId << "exists!\n";
			continue;
		}

		// align the read to the somatic contig
		splitReadAlignment bp1ContigSR;
		if (svAlignInfo.bp1ContigReversed)
			bp1ContigSR.align(readSeqRC, svAlignInfo.contigSeq, svAlignInfo.bp1ContigOffset);
		else
			bp1ContigSR.align(readSeq, svAlignInfo.contigSeq, svAlignInfo.bp1ContigOffset);
		splitReadAlignment bp2ContigSR;
		if (svAlignInfo.bp2ContigReversed)
			bp2ContigSR.align(readSeqRC, svAlignInfo.contigSeq, svAlignInfo.bp2ContigOffset);
		else
			bp2ContigSR.align(readSeq, svAlignInfo.contigSeq, svAlignInfo.bp2ContigOffset);

		// align the read to reference regions
		splitReadAlignment bp1RefSR;
		bp1RefSR.align(readSeq, svAlignInfo.bp1RefSeq, svAlignInfo.bp1RefOffset);
		splitReadAlignment bp2RefSR;
		bp2RefSR.align(readSeq, svAlignInfo.bp2RefSeq, svAlignInfo.bp2RefOffset);

		// scoring
		float bp1ContigEvidence(0);
		float bp2ContigEvidence(0);
		float bp1RefEvidence(0);
		float bp2RefEvidence(0);
		if (bp1ContigSR.has_evidence()) bp1ContigEvidence = bp1ContigSR.get_evidence();
		if (bp2ContigSR.has_evidence()) bp2ContigEvidence = bp2ContigSR.get_evidence();
		if (bp1RefSR.has_evidence()) bp1RefEvidence = bp1RefSR.get_evidence();
		if (bp2RefSR.has_evidence()) bp2RefEvidence = bp2RefSR.get_evidence();

		const float contigEvidence = std::max(bp1ContigEvidence, bp2ContigEvidence);
		const float refEvidence = std::max(bp1RefEvidence, bp2RefEvidence);

		if ((bp1ContigSR.has_evidence()) || (bp1ContigSR.has_evidence()))
		{
			sample.contigSRCount++;
			sample.contigSREvidence += contigEvidence;
			sample.contigSRMapQ += readMapQ * readMapQ;

#ifdef DEBUG_SVS
			/*
			log_os << "bp1Contig \n";
			log_os << bp1ContigSR << "\n";
			log_os << "bp2Contig \n";
			log_os << bp2ContigSR << "\n";
			*/

#endif
		}
		if ((bp1RefSR.has_evidence()) || (bp2RefSR.has_evidence()))
		{
			sample.refSRCount++;
			sample.refSREvidence += refEvidence;
			sample.refSRMapQ += readMapQ * readMapQ;

#ifdef DEBUG_SVS
			/*
			log_os << "bp1Ref \n";
			log_os << bp1RefSR << "\n";
			log_os << "bp2Ref \n";
			log_os << bp2RefSR << "\n";
			*/
#endif
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

    // extract SV alignment info for split read evidence
    const SVAlignmentInfo SVAlignInfo(sv, assemblyData);

    // first exercise -- just count the sample assignment of the pairs we already have:
    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const bool isTumor(_isAlignmentTumor[bamIndex]);
        SVSampleInfo& sample(isTumor ? ssInfo.tumor : ssInfo.normal);

        // consider 2-locus events first
        // TODO: to add local assembly later
        // apply the split-read scoring, only when the values of max depth are reasonable;
        // otherwise, the read map may blow out
        if ((assemblyData.isSpanning) &&
        	(ssInfo.bp1MaxDepth <= 2*_dFilter.maxDepth(sv.bp1.interval.tid)) &&
        	(ssInfo.bp2MaxDepth <= 2*_dFilter.maxDepth(sv.bp2.interval.tid)))
        {
        	read_map_t readMap;

        	streamPtr& bamPtr(_bamStreams[bamIndex]);
			bam_streamer& read_stream(*bamPtr);
			// scoring split reads overlapping bp1
			scoreSplitReads(true, sv.bp1, SVAlignInfo, readMap,
					        read_stream, sample);
			// scoring split reads overlapping bp2
			scoreSplitReads(false, sv.bp2, SVAlignInfo, readMap,
			        		read_stream, sample);

#ifdef DEBUG_SVS
			log_os << "\nbam is tumor = " << isTumor << "\n";
			log_os << "tumor contig SP count: " << ssInfo.tumor.contigSRCount << "\n";
			log_os << "tumor contig SP evidence: " << ssInfo.tumor.contigSREvidence << "\n";
			log_os << "normal contig SP count: " << ssInfo.normal.contigSRCount << "\n";
			log_os << "normal contig SP evidence: " << ssInfo.normal.contigSREvidence << "\n";

			log_os << "tumor ref SP count: " << ssInfo.tumor.refSRCount << "\n";
			log_os << "tumor ref SP evidence: " << ssInfo.tumor.refSREvidence << "\n";
			log_os << "normal ref SP count: " << ssInfo.normal.refSRCount << "\n";
			log_os << "normal ref SP evidence: " << ssInfo.normal.refSREvidence << "\n";
#endif
        }


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

    // root mean square
    if (ssInfo.tumor.contigSRCount > 0)
    	ssInfo.tumor.contigSRMapQ = sqrt(ssInfo.tumor.contigSRMapQ / (float)ssInfo.tumor.contigSRCount);
    if (ssInfo.tumor.refSRCount > 0)
    	ssInfo.tumor.refSRMapQ = sqrt(ssInfo.tumor.refSRMapQ / (float)ssInfo.tumor.refSRCount);
    if (ssInfo.normal.contigSRCount > 0)
    	ssInfo.normal.contigSRMapQ = sqrt(ssInfo.normal.contigSRMapQ / (float)ssInfo.normal.contigSRCount);
    if (ssInfo.normal.refSRCount > 0)
    	ssInfo.normal.refSRMapQ = sqrt(ssInfo.normal.refSRMapQ / (float)ssInfo.normal.refSRCount);

#ifdef DEBUG_SVS
    log_os << "\nfinally...\n";
    log_os << "tumor contig SP count: " << ssInfo.tumor.contigSRCount << "\n";
    log_os << "tumor contig SP evidence: " << ssInfo.tumor.contigSREvidence << "\n";
    log_os << "tumor contig SP_mapQ: " << ssInfo.tumor.contigSRMapQ << "\n";
    log_os << "normal contig SP count: " << ssInfo.normal.contigSRCount << "\n";
    log_os << "normal contig SP evidence: " << ssInfo.normal.contigSREvidence << "\n";
    log_os << "normal contig SP_mapQ: " << ssInfo.normal.contigSRMapQ << "\n";

    log_os << "tumor ref SP count: " << ssInfo.tumor.refSRCount << "\n";
    log_os << "tumor ref SP evidence: " << ssInfo.tumor.refSREvidence << "\n";
    log_os << "tumor ref SP_mapQ: " << ssInfo.tumor.refSRMapQ << "\n";
    log_os << "normal ref SP count: " << ssInfo.normal.refSRCount << "\n";
    log_os << "normal ref SP evidence: " << ssInfo.normal.refSREvidence << "\n";
    log_os << "normal ref SP_mapQ: " << ssInfo.normal.refSRMapQ << "\n";
#endif



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




