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


#include "ReadGroupStats.hh"

#include "blt_util/bam_streamer.hh"
#include "blt_util/log.hh"

#include "boost/foreach.hpp"

#include <vector>
#include <iostream>

//#define DEBUG_RPS



static
bool
isStatSetMatch(const SizeDistribution& pss1,
               const SizeDistribution& pss2)
{
    static const float cdfPrecision(0.005);

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
getRelOrient(const bam_record& br)
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
    static const unsigned statsCheckCnt(100000);
    static const unsigned maxPosCount(1);

    bam_streamer read_stream(statsBamFile.c_str());

    const bam_header_t& header(* read_stream.get_header());
    const int32_t nChrom(header.n_targets);
    std::vector<int32_t> chromSize(nChrom,0);
    std::vector<int32_t> chromHighestPos(nChrom,-1);
    for (int32_t i(0); i<nChrom; ++i)
    {
        chromSize[i] = (header.target_len[i]);
    }

    bool isConverged(false);
    bool isStopEstimation(false);
    bool isFirstEstimation(true);
    SizeDistribution oldFragSize;

    unsigned recordCnts(0);
    unsigned posCount(0);
    bool isPairTypeSet(false);
    bool isActiveChrom(true);

    while (isActiveChrom && (!isStopEstimation))
    {
        isActiveChrom=false;
        for (int32_t i(0); i<nChrom; ++i)
        {
            if (isStopEstimation) break;

            const int32_t startPos(chromHighestPos[i]+1);
#ifdef DEBUG_RPS
            std::cerr << "INFO: Stats requesting bam region starting from: chrid: " << i << " start: " << startPos << "\n";
#endif
            read_stream.set_new_region(i,startPos,chromSize[i]);
            while (read_stream.next())
            {
                const bam_record& al(*(read_stream.get_record_ptr()));
                if (al.pos()<startPos) continue;

                if (al.pos()!=chromHighestPos[i])
                {
                    posCount=0;
                }

                chromHighestPos[i]=al.pos();
                isActiveChrom=true;

                if (! (al.is_paired() && al.is_proper_pair())) continue;
                if (al.map_qual()==0) continue;

                // sample each read pair once by sampling stats from
                // upstream read only:
                if (al.pos()<al.mate_pos()) continue;

                // to prevent high-depth pileups from overly biasing the
                // read stats, we only take maxPosCount read pairs from each start
                // pos:
                if (posCount>=maxPosCount) continue;
                posCount++;

                ++recordCnts;

                // Assuming only two reads per fragment - based on bamtools.
                const unsigned readNum(al.is_first() ? 1 : 2);
                assert(al.is_second() == (readNum == 2));

                if (! isPairTypeSet)
                {
                    // TODO: does orientation need to be averaged over several observations?
                    relOrients = getRelOrient(al);
                    isPairTypeSet=true;
                }

                const unsigned currFragSize(std::abs(al.template_size()));
                fragStats.addObservation(currFragSize);

                if ((recordCnts % statsCheckCnt) != 0) continue;

#ifdef DEBUG_RPS
                log_os << "INFO: Checking stats convergence at record count : " << recordCnts << "'\n"
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

                if (isConverged || (recordCnts>5000000)) isStopEstimation=true;

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
        else if ((recordCnts % statsCheckCnt) != 0)
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
}
