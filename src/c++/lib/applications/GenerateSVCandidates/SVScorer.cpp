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
#include "blt_util/bam_record_util.hh"
#include "blt_util/log.hh"
#include "blt_util/ReadKey.hh"
#include "common/Exceptions.hh"
#include "manta/ReadGroupStatsSet.hh"
#include "boost/foreach.hpp"

#include <algorithm>
#include <iostream>
#include <string>

//#define DEBUG_SVS


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



typedef std::set<ReadKey> read_map_t;



static
void
scoreSplitReads(
    const bool isBp1,
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
        const ReadKey readKey(bamRead);

        // map all reads for bp1
        // then for bp2, skip reads that have been considered for bp1 to avoid double-counting
        if (isBp1)
        {
            readMap.insert(readKey);
        }
        if (readMap.find(readKey) != readMap.end())
        {
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

        if ((bp1ContigSR.has_evidence()) || (bp2ContigSR.has_evidence()))
        {
            sample.contigSRCount++;
            sample.contigSREvidence += contigEvidence;
            sample.contigSRMapQ += readMapQ * readMapQ;

#ifdef DEBUG_SVS
            /*
            log_os << "bp1Contig \n";
            log_os << bp1ContigSR;
            log_os << "bp2Contig \n";
            log_os << bp2ContigSR;
            log_os << "contigEvidence = " << contigEvidence << "\n";
            log_os << "accumulated contigEvidence = " << sample.contigSREvidence << "\n";
            log_os << "contigCount = " << sample.contigSRCount << "\n\n";
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
            log_os << bp1RefSR;
            log_os << "bp2Ref \n";
            log_os << bp2RefSR;
            log_os << "refEvidence = " << refEvidence << "\n";
            log_os << "accumulated refEvidence = " << sample.refSREvidence << "\n";
            log_os << "refCount = " << sample.refSRCount << "\n\n";
            */
#endif
        }
    }
}



/// get reference allele support at a single breakend:
///
void
SVScorer::
getSVRefPairSupport(
    const SVBreakend& bp,
    SomaticSVScoreInfo& ssInfo,
    const bool isBp1)
{
    /// search for all read pairs supporting the reference allele
    ///
    /// APPROXIMATION: for imprecise and precise variants treat the breakend locations as the center of the
    ///  breakend interval.
    ///
    /// TODO: improve on the approx above
    ///

    /// TODO: track read key to account for overlap of spanning and split read evidence

    const pos_t centerPos(bp.interval.range.center_pos());

    
    /// we're interested in any fragments which cross center pos with at least N bases of support on each side
    /// (note this definition is certain to overlap the split read definition whenever N is less than the read length
    
    static const pos_t minFragSupport(80);
    
    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const bool isTumor(_isAlignmentTumor[bamIndex]);
        SVSampleInfo& sample(isTumor ? ssInfo.tumor : ssInfo.normal);

        bam_streamer& bamStream(*_bamStreams[bamIndex]);

        /// set the search range around centerPos so that we can get any fragments at the Xth percentile length or smaller which could have
        /// min Fragsupport
        const SVLocusScanner::Range& pRange(_readScanner.getEvidencePairRange(bamIndex));
        const unsigned minFrag(pRange.min);
        const unsigned maxFrag(pRange.max);

        const uint32_t beginPos(centerPos-maxFrag+minFragSupport);
        const uint32_t endPos(centerPos+maxFrag-minFragSupport+1);
        
        if (beginPos >= endPos) continue;

        // set bam stream to new search interval:
        bamStream.set_new_region(bp.interval.tid, beginPos, endPos);

        while (bamStream.next())
        {
            const bam_record& bamRead(*(bamStream.get_record_ptr()));

            if (_readScanner.isReadFiltered(bamRead)) continue;
            if (bamRead.is_unmapped() || bamRead.is_mate_unmapped()) continue;

            /// check for standard innie orientation:
            if (! is_innie_pair(bamRead)) continue;

            /// check if fragment is too big or too small:
            const unsigned tSize(std::abs(bamRead.template_size()));
            if (tSize < minFrag) continue;
            if (tSize > maxFrag) continue;

            // count only from the down stream read unless the mate-pos goes past center-pos
            const bool isLeftMost(bamRead.pos() < bamRead.mate_pos());
            const bool isMateBeforeCenter(bamRead.mate_pos() < centerPos);

            if ( isLeftMost && isMateBeforeCenter ) continue;
            if ( (!isLeftMost) && (!isMateBeforeCenter) ) continue;

            // get fragment range:
            pos_t fragBegin(0);
            pos_t fragEnd(0);
            if (isLeftMost)
            {
                fragBegin=bamRead.pos()-1;
                fragEnd=fragBegin+bamRead.template_size();
            }
            else
            {
                fragBegin=bamRead.mate_pos()-1;
                fragEnd=fragBegin-bamRead.template_size();
            }

            if (fragBegin > fragEnd)
            {
                using namespace illumina::common;

                std::ostringstream oss;
                oss << "ERROR: Failed to parse fragment range from bam record. Frag begin,end: " << fragBegin << " " << fragEnd << " bamRecord: " << bamRead << "\n";
                BOOST_THROW_EXCEPTION(LogicException(oss.str()));
            }

            const pos_t fragOverlap(std::min((centerPos-fragBegin), (fragEnd-centerPos)));

            if (fragOverlap < minFragSupport) continue;

            if (isBp1)
            {
                sample.refAlleleBp1SpanPairs++;
            }
            else
            {
                sample.refAlleleBp2SpanPairs++;
            }
        }
    }
}



// make final interpretation of reference support as the minimum breakend support:
static
void
finishSamplePairSupport(
    SVSampleInfo& sample)
{
    sample.refAlleleSpanPairs = std::min(sample.refAlleleBp1SpanPairs, sample.refAlleleBp2SpanPairs);
}



void
SVScorer::
getSVRefPairSupport(
    const SVCandidate& sv,
    SomaticSVScoreInfo& ssInfo)
{    
    getSVRefPairSupport(sv.bp1, ssInfo, true);
    getSVRefPairSupport(sv.bp2, ssInfo, false);

    finishSamplePairSupport(ssInfo.tumor);
    finishSamplePairSupport(ssInfo.normal);
}



/// return rms given sum of squares
static
float
finishRms(
    const float sumSqr,
    const unsigned count)
{
    if (count == 0) return 0.;
    return std::sqrt(sumSqr / static_cast<float>(count));
}



/// make final split read computations after bam scanning is finished:
static
void
finishSampleSRData(
    SVSampleInfo& sample)
{
    // root mean square
    sample.contigSRMapQ = finishRms(sample.contigSRMapQ, sample.contigSRCount);
    sample.refSRMapQ = finishRms(sample.refSRMapQ, sample.refSRCount);
}



void
SVScorer::
getSVSplitReadSupport(
    const SVCandidateAssemblyData& assemblyData,
    const SVCandidate& sv,
    SomaticSVScoreInfo& ssInfo)
{
    static const unsigned maxDepthSRFactor(2); ///< at what multiple of the maxDepth do we skip split read analysis?
    const bool isSkipSRSearchDepth(
        (ssInfo.bp1MaxDepth > maxDepthSRFactor*_dFilter.maxDepth(sv.bp1.interval.tid)) ||
        (ssInfo.bp2MaxDepth > maxDepthSRFactor*_dFilter.maxDepth(sv.bp2.interval.tid)));

    // apply the split-read scoring, only when:
    // 1) the SV is precise, i.e. has successful somatic contigs;
    // 2) the values of max depth are reasonable (otherwise, the read map may blow out).
    //
    // consider 2-locus events first
    // TODO: to add local assembly later
    //
    const bool isSkipSRSearch(
            (! assemblyData.isSpanning) ||
            (sv.isImprecise()) ||
            (isSkipSRSearchDepth));

    if(isSkipSRSearch) return;

    // Get Data on standard read pairs crossing the two breakends,

    // extract SV alignment info for split read evidence
    const SVAlignmentInfo SVAlignInfo(sv, assemblyData);


    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const bool isTumor(_isAlignmentTumor[bamIndex]);
        SVSampleInfo& sample(isTumor ? ssInfo.tumor : ssInfo.normal);
        bam_streamer& bamStream(*_bamStreams[bamIndex]);

        read_map_t readMap;

        // scoring split reads overlapping bp1
        scoreSplitReads(true, sv.bp1, SVAlignInfo, readMap,
                        bamStream, sample);
        // scoring split reads overlapping bp2
        scoreSplitReads(false, sv.bp2, SVAlignInfo, readMap,
                        bamStream, sample);
    }

    finishSampleSRData(ssInfo.tumor);
    finishSampleSRData(ssInfo.normal);

#ifdef DEBUG_SVS
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

    // count the split reads supporting the ref and alt alleles in each sample
    //
    getSVSplitReadSupport(assemblyData, sv, ssInfo);

#ifdef DEBUG_SVS
    log_os << SVAlignInfo << "\n";
#endif

    // count the read pairs supporting the alternate allele in each sample, using data we already produced during candidate generation:
    //
    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const bool isTumor(_isAlignmentTumor[bamIndex]);
        SVSampleInfo& sample(isTumor ? ssInfo.tumor : ssInfo.normal);

        const SVCandidateSetReadPairSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
        BOOST_FOREACH(const SVCandidateSetReadPair& pair, svDataGroup)
        {
            if (0 == std::count(pair.svIndex.begin(),pair.svIndex.end(), sv.candidateIndex)) continue;

            if (pair.read1.isSet())
            {
                sample.altAlleleBp1SpanReads += 1;
            }
            if (pair.read2.isSet())
            {
                sample.altAlleleBp2SpanReads += 1;
            }
            if (pair.read1.isSet() && pair.read2.isSet())
            {
                sample.altAlleleSpanPairs += 1;
            }
        }
    }

    // count the read pairs supporting the reference allele on each breakend in each sample:
    //
    getSVRefPairSupport(sv,ssInfo);


    //
    // apply filters
    //
    if (_dFilter.isMaxDepthFilter())
    {
        // apply maxdepth filter if either of the breakpoints exceeds the maximum depth:
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
    if (ssInfo.normal.altAlleleSpanPairs > 1) isSomatic=false;

    if (isSomatic)
    {
        const bool lowPairSupport(ssInfo.tumor.altAlleleSpanPairs < 6);
        const bool lowSingleSupport((ssInfo.tumor.altAlleleBp1SpanReads < 14) || (ssInfo.tumor.altAlleleBp2SpanReads < 14));
        const bool highSingleContam((ssInfo.normal.altAlleleBp1SpanReads > 1) || (ssInfo.normal.altAlleleBp2SpanReads > 1));

        /// allow single pair support to rescue an SV only if the evidence looks REALLY good:
        if (lowPairSupport && (lowSingleSupport || highSingleContam))
            isSomatic=false;
    }

    if (isSomatic)
    {
        if (ssInfo.normal.altAlleleSpanPairs)
        {
            const double ratio(static_cast<double>(ssInfo.tumor.altAlleleSpanPairs)/static_cast<double>(ssInfo.normal.altAlleleSpanPairs));
            if (ratio<9)
            {
                isSomatic=false;
            }
        }
        if (ssInfo.normal.altAlleleBp1SpanReads)
        {
            const double ratio(static_cast<double>(ssInfo.tumor.altAlleleBp1SpanReads)/static_cast<double>(ssInfo.normal.altAlleleBp1SpanReads));
            if (ratio<9)
            {
                isSomatic=false;
            }
        }
        if (ssInfo.normal.altAlleleBp2SpanReads)
        {
            const double ratio(static_cast<double>(ssInfo.tumor.altAlleleBp2SpanReads)/static_cast<double>(ssInfo.normal.altAlleleBp2SpanReads));
            if (ratio<9)
            {
                isSomatic=false;
            }
        }
    }

    if (isSomatic) ssInfo.somaticScore=60;
}




