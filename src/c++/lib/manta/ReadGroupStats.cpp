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
/// \author Bret Barnes, Xiaoyu Chen
///

#include "ReadGroupStats.hh"

#include "blt_util/align_path_bam_util.hh"
#include "blt_util/bam_record_util.hh"
#include "blt_util/bam_streamer.hh"
#include "blt_util/log.hh"
#include "blt_util/ReadKey.hh"
#include "manta/SVLocusScanner.hh"

#include "boost/foreach.hpp"

#include <iostream>
#include <set>
#include <vector>

//#define DEBUG_RPS



static
bool
isStatSetMatch(
    const SizeDistribution& pss1,
    const SizeDistribution& pss2)
{
    static const float cdfPrecision(0.001);

    for (float prob(0.05); prob < 1; prob += 0.1)
    {
        // check if percentile values equal
        if (std::abs(pss1.quantile(prob) - pss2.quantile(prob)) >= 1)
        {
            return false;
        }

        // check the convergence of fragsize cdf
        const int fragSize(pss2.quantile(prob));
        if (std::abs(pss1.cdf(fragSize) - pss2.cdf(fragSize)) >= cdfPrecision)
        {
            return false;
        }
    }

    return true;
}



/// This produces a useful result only when both reads align to the same
/// chromosome.
static
ReadPairOrient
getRelOrient(
    const bam_record& br)
{
    pos_t pos1 = br.pos();
    bool is_fwd_strand1 = br.is_fwd_strand();
    pos_t pos2 = br.mate_pos();
    bool is_fwd_strand2 = br.is_mate_fwd_strand();

    if (! br.is_first())
    {
        std::swap(pos1,pos2);
        std::swap(is_fwd_strand1,is_fwd_strand2);
    }

    ReadPairOrient rpo;
    rpo.setVal(PAIR_ORIENT::get_index(pos1,is_fwd_strand1,pos2,is_fwd_strand2));
    return rpo;
}



/* ----- ----- ----- ----- ----- -----
 *
 * ----- ReadPairStats  -----
 *
 * ----- ----- ----- ----- ----- ----- */

// set read pair statistics from a bam reader object:
//
ReadGroupStats::
ReadGroupStats(const std::string& statsBamFile)
{
    typedef std::set<ReadKey> mateMap_t;

    static const unsigned statsCheckCnt(100000);
    static const unsigned maxPosCount(1);

    bam_streamer read_stream(statsBamFile.c_str());

    const bam_header_t& header(* read_stream.get_header());
    const int32_t chromCount(header.n_targets);
    std::vector<int32_t> chromSize(chromCount,0);
    std::vector<int32_t> chromHighestPos(chromCount,-1);
    for (int32_t i(0); i<chromCount; ++i)
    {
        chromSize[i] = (header.target_len[i]);
    }

    // cache-variables -- these do not need to be stored across loop iterations, we just save sys calls by keeping them here
    ALIGNPATH::path_t apath;
    mateMap_t goodMates;

    bool isConverged(false);
    bool isStopEstimation(false);
    bool isFirstEstimation(true);
    SizeDistribution oldFragSize;

    unsigned recordCount(0);
    unsigned posCount(0);
    bool isPairTypeSet(false);
    bool isActiveChrom(true);

    while (isActiveChrom && (!isStopEstimation))
    {
        isActiveChrom=false;
        for (int32_t chromIndex(0); chromIndex<chromCount; ++chromIndex)
        {
            if (isStopEstimation) break;

            goodMates.clear();

            const int32_t startPos(chromHighestPos[chromIndex]+1);
#ifdef DEBUG_RPS
            std::cerr << "INFO: Stats requesting bam region starting from: chrid: " << chromIndex << " start: " << startPos << "\n";
#endif
            read_stream.set_new_region(chromIndex,startPos,chromSize[chromIndex]);
            while (read_stream.next())
            {
                const bam_record& bamRead(*(read_stream.get_record_ptr()));
                if (bamRead.pos()<startPos) continue;

                if (bamRead.pos() != chromHighestPos[chromIndex])
                {
                    posCount=0;
                }

                chromHighestPos[chromIndex]=bamRead.pos();
                isActiveChrom=true;

                // filter common categories of undesirable reads:
                if (SVLocusScanner::isReadFilteredCore(bamRead)) continue;

                // filter mapped innies on the same chrom
                //
                // note we don't rely on the proper pair bit because this already contains an arbitrary length filter AND subjects the method to
                // aligner specific variation
                //
                // TODO: ..note this locks-in standard ilmn orientation -- okay for now but this function needs major re-arrangement for mate-pair support,
                // we could still keep independence from each aligner's proper pair decisions by estimating a fragment distro for each orientation
                // and only keeping the one with the most samples
                if (! is_innie_pair(bamRead)) continue;
                if (bamRead.map_qual()==0) continue;

                // filter any split reads with an SA tag:
                static const char SAtag[] = {'S','A'};
                if (NULL != bamRead.get_string_tag(SAtag)) continue;

                bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), apath);

                {
                    // use only the most conservative alignments to generate fragment stats --
                    // filter reads containing any cigar types besides MATCH:
                    bool isBadAlign(false);
                    BOOST_FOREACH(const ALIGNPATH::path_segment& ps, apath)
                    {
                        if (! ALIGNPATH::is_segment_align_match(ps.type))
                        {
                            isBadAlign=true;
                            break;
                        }
                    }
                    if (isBadAlign) continue;
                }

                // sample each read pair once by sampling stats from
                // downstream read only, or whichever read is encountered
                // second if they're at the same position:
                bool isDownstream(bamRead.pos() > bamRead.mate_pos());

                if (bamRead.pos() == bamRead.mate_pos())
                {
                    /// in this case we consider the first read we run into to be "downstream",
                    /// the first read is determined by looking in the mate set:
                    const int mateReadNo( bamRead.is_first() ? 2 : 1);
                    const ReadKey mateKey(bamRead.qname(), mateReadNo);

                    mateMap_t::iterator i(goodMates.find(mateKey));

                    if (i != goodMates.end()) isDownstream = true;
                }

                if (! isDownstream)
                {
                    // to prevent high-depth pileups from overly biasing the
                    // read stats, we only take maxPosCount read pairs from each start
                    // pos:
                    if (posCount>=maxPosCount) continue;

                    /// crude mechanism to manage total set memory
                    static const unsigned maxMateSetSize(100000);
                    if (goodMates.size() > maxMateSetSize) goodMates.clear();

                    const ReadKey key(bamRead);
                    goodMates.insert(key);

                    ++posCount;

                    continue;
                }
                else
                {
                    const int mateReadNo( bamRead.is_first() ? 2 : 1);
                    const ReadKey mateKey(bamRead.qname(), mateReadNo);

                    mateMap_t::iterator i(goodMates.find(mateKey));

                    if (i == goodMates.end()) continue;

                    goodMates.erase(i);
                }

                // made it through all filters!
                ++recordCount;

                // Assuming only two reads per fragment - based on bamtools.
                const unsigned readNum(bamRead.is_first() ? 1 : 2);
                assert(bamRead.is_second() == (readNum == 2));

                if (bamRead.pos() != bamRead.mate_pos())
                {
                    if (! isPairTypeSet)
                    {
                        // TODO: does orientation need to be averaged over several observations?
                        relOrients = getRelOrient(bamRead);
                        isPairTypeSet=true;
                    }
                }

                unsigned currFragSize(std::abs(bamRead.template_size()));

                // reduce fragsize resolution for very large sizes:
                // (large sizes are uncommon -- this doesn't need to be clever/fast)
                {
                    unsigned steps(0);
                    while (currFragSize>1000)
                    {
                        currFragSize /= 10;
                        steps++;
                    }
                    for (unsigned stepIndex(0); stepIndex<steps; ++stepIndex) currFragSize *= 10;
                }

                fragStats.addObservation(currFragSize);

                if ((recordCount % statsCheckCnt) != 0) continue;

#ifdef DEBUG_RPS
                log_os << "INFO: Checking stats convergence at record count : " << recordCount << "'\n"
                       << "INFO: Stats before convergence check: ";
                //write(log_os);
                log_os << "\n";
#endif

                // check convergence
                if (isFirstEstimation)
                {
                    isFirstEstimation = false;
                }
                else
                {
                    isConverged=isStatSetMatch(oldFragSize, fragStats);
                }

                oldFragSize = fragStats;

                static const unsigned maxRecordCount(5000000);
                if (isConverged || (recordCount>maxRecordCount)) isStopEstimation=true;

                // break from reading the current chromosome
                break;
            }
        }
    }

    if (!isConverged)
    {
        if (fragStats.totalObservations() <1000)
        {
            log_os << "ERROR: Can't generate pair statistics for BAM file " << statsBamFile << "\n";
            log_os << "\tTotal observed read pairs: " << fragStats.totalObservations() << "\n";
            exit(EXIT_FAILURE);
        }
        else if ((recordCount % statsCheckCnt) != 0)
        {
            if (! isFirstEstimation)
            {
                isConverged=isStatSetMatch(oldFragSize, fragStats);
            }
        }

        if (!isConverged)
        {
            log_os << "WARNING: read pair statistics did not converge\n";
        }
    }

    // final step before saving is to cut-off the extreme end of the fragment size distribution, this
    // is similar the some aligners proper-pair bit definition of (3x the standard mean, etc.)
    static const float filterQuant(0.9995);
    fragStats.filterObservationsOverQuantile(filterQuant);
}
