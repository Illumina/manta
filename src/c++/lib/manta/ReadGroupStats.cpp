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
#include "blt_util/parse_util.hh"

#include <boost/foreach.hpp>

#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>

#define DEBUG_RPS

/* ----- ----- ----- ----- ----- -----
 *
 * ----- Auxiliary functions     -----
 *
 * ----- ----- ----- ----- ----- ----- */
static
void
writeFragSizeHashItem(std::ostream& os, PairStatSet::hash_map_fragment fragmentSizeHash, int k)
{
    os << k << " -> <"
       << fragmentSizeHash.find(k)->second.first
       << ", "
       << fragmentSizeHash.find(k)->second.second
       << ">\n";
}


static
void
populateCdfQuantiles(PairStatSet::hash_map_fragment& fragmentSizeHash,
                     std::vector<int> fragmentSizes,
                     int numOfFragSize, int totalCount,
                     int quantileNum, float* quantiles)
{
    int fillBase = 0;
    float cumulative = 0;
    for (int s=0; s<numOfFragSize; s++)
    {
        int fs = fragmentSizes[s];
        int count = fragmentSizeHash.find(fs)->second.first;
        float freq = count / (float)totalCount;

        cumulative += freq;
        // update the hash map with cdf
        fragmentSizeHash.find(fs)->second.second = cumulative;

        int fillNext = rint(cumulative * quantileNum);
        for (int q = fillBase; q < fillNext; q++)
            quantiles[q] = fs;
        fillBase = fillNext;
    }
}


static
bool
isStatSetMatch(const PairStatSet& pss1,
               const PairStatSet& pss2)
{
    static const double statsPrecision(0.005);

    float p = 0.05;
    float delta = 0.1;
    while (p < 1)
    {
        // check if percentile values equal
        int b = p * pss2.quantileNum - 1;
        if (std::abs(pss1.quantiles[b] - pss2.quantiles[b])>=1)
            return false;

        // check the convergence of cdf of the median value
        int medFragSize = pss2.quantiles[b];
        if (std::abs(pss1.cdf(medFragSize) - pss2.cdf(medFragSize)) >= statsPrecision)
            return false;

        p += delta;
    }

    return true;
}

static
int
binarySearch(int vectorLen, std::vector<int> sortedVec, int value)
{
    int ret = -1;

    int lowIx = 0;
    int highIx = vectorLen;
    while (lowIx + 1 < highIx)
    {
        int midIx = (highIx + lowIx) / 2;
        if (sortedVec[midIx] > value)
            highIx = midIx;
        else
            lowIx = midIx;
    }

    if (value >= sortedVec[lowIx])
        ret = sortedVec[lowIx];

    return ret;
}



// This produces a useful result only when both reads align to the same
// chromosome.
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
 * ----- PairStatSet  -----
 * ----- ----- ----- ----- ----- ----- */
bool
PairStatSet::
calcStats()
{
#ifdef DEBUG_RPS
    std::cerr<<"Calculating stats...\n";
#endif
    // calculate statistics from hashed insert sizes
    numOfFragSize = fragmentSizeHash.size();
#ifdef DEBUG_RPS
    std::cerr<<"numOfFragSize="<<numOfFragSize<<"\n";
#endif
    if (numOfFragSize == 0)
        return false;

    // clean the vector
    fragmentSizes.clear();
    // populate the vector of fragment sizes
    BOOST_FOREACH (hash_map_fragment::value_type& hashItem, fragmentSizeHash)
    fragmentSizes.push_back(hashItem.first);
    // sort all insert sizes
    std::sort(fragmentSizes.begin(), fragmentSizes.end());

    // populate the array of quantiles
    populateCdfQuantiles(fragmentSizeHash, fragmentSizes, numOfFragSize,
                         totalCount, quantileNum, quantiles);


    return true;
}



int
PairStatSet::
quantile(const float p) const
{
    int insertSize = 0;

    int bin = rint(p * quantileNum) - 1;
    if ((bin >= 0) && (bin < quantileNum))
        insertSize = quantiles[bin];

    return insertSize;
}

float
PairStatSet::
cdf(const int fs) const
{
    float cumProb = 0;

    if (fragmentSizeHash.find(fs) != fragmentSizeHash.end())
        cumProb = fragmentSizeHash.find(fs)->second.second;
    else
    {
        int estimated = binarySearch(numOfFragSize, fragmentSizes, fs);
        if (estimated > -1)
            cumProb = fragmentSizeHash.find(fs)->second.second;
    }

    return cumProb;
}

std::ostream&
operator<<(std::ostream& os, const PairStatSet& pss)
{
    os << pss.totalCount << '\t'
       << pss.numOfFragSize << '\n';

#ifdef DEBUG_RPS
    // testing...
    os << "cdf(0)=" << pss.cdf(0)<<"\n"
       << "cdf(100)=" << pss.cdf(100)<<"\n"
       << "cdf(200)=" << pss.cdf(200)<<"\n"
       << "cdf(300)=" << pss.cdf(300)<<"\n"
       << "cdf(400)=" << pss.cdf(400)<<"\n"
       << "cdf(500)=" << pss.cdf(500)<<"\n"
       << "cdf(600)=" << pss.cdf(600)<<"\n"
       << "cdf(700)=" << pss.cdf(700)<<"\n"
       << "cdf(800)=" << pss.cdf(800)<<"\n"
       << "cdf(900)=" << pss.cdf(900)<<"\n"
       << "cdf(1000)=" << pss.cdf(1000)<<"\n";

    os << "quantile(0.12)=" << pss.quantile(0.12)<<"\n"
       << "quantile(0.1201)=" << pss.quantile(0.1201)<<"\n"
       << "quantile(0.1205)=" << pss.quantile(0.1205)<<"\n"
       << "quantile(0.5)=" << pss.quantile(0.5)<<"\n"
       << "quantile(0.5001)=" << pss.quantile(0.5001)<<"\n"
       << "quantile(0.5006)=" << pss.quantile(0.5006)<<"\n"
       << "quantile(0.99)=" << pss.quantile(0.99)<<"\n"
       << "quantile(0.9991)=" << pss.quantile(0.9991)<<"\n"
       << "quantile(0.9998)=" << pss.quantile(0.9998)<<"\n";
#endif

    return os;
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
    PairStatSet oldFragSize;

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
                const unsigned int readNum(al.is_first() ? 1 : 2);
                assert(al.is_second() == (readNum == 2));

                if (! isPairTypeSet)
                {
                    // TODO: does orientation need to be averaged over several observations?
                    relOrients = getRelOrient(al);
                    isPairTypeSet=true;
                }

                fragSize.totalCount++;
                int currFragSize = std::abs(al.template_size());
                if (fragSize.fragmentSizeHash.find(currFragSize) == fragSize.fragmentSizeHash.end())
                    // initialize the count
                    fragSize.fragmentSizeHash[currFragSize] = std::make_pair(1, 0);
                else
                    // increase the count
                    fragSize.fragmentSizeHash[currFragSize].first++;

                if ((recordCnts % statsCheckCnt) != 0) continue;

#ifdef DEBUG_RPS
                log_os << "INFO: Checking stats convergence at record count : " << recordCnts << "'\n"
                       << "INFO: Stats before convergence check: ";
                //write(log_os);
                log_os << "\n";
#endif

                fragSize.calcStats();
                // check convergence
                if (isFirstEstimation)
                    isFirstEstimation = false;
                else
                    isConverged=isStatSetMatch(oldFragSize, fragSize);

                oldFragSize = fragSize;

                if (isConverged || (recordCnts>5000000))
                    isStopEstimation=true;

                // break from reading the current chromosome
                break;
            }
        }
    }

    if (!isConverged)
    {
        if (fragSize.totalCount <1000)
        {
            log_os << "ERROR: Can't generate pair statistics for BAM file " << statsBamFile << "\n";
            log_os << "\tTotal observed read pairs: " << fragSize.totalCount << "\n";
            exit(EXIT_FAILURE);
        }
        else if ((recordCnts % statsCheckCnt) != 0)
        {
            // caculate for
            fragSize.calcStats();
            // check convergence
            if (isFirstEstimation)
                isFirstEstimation = false;
            else
                isConverged=isStatSetMatch(oldFragSize, fragSize);
        }

        if (!isConverged)
            log_os << "WARNING: read pair statistics did not converge\n";
    }
}


void
ReadGroupStats::
write(std::ostream& os) const
{
    os << relOrients << "\t"
       << fragSize;
}

