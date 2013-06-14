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

#include <cmath>

#include <iostream>
#include <vector>

//#define DEBUG_RPS


// Stats file data format
const unsigned HEAD_FILL_IDX = 0;
const unsigned HEAD_SRC_IDX  = 1;
const unsigned HEAD_NAME_IDX = 2;

// Stats file data format
const unsigned STAT_SOURCE_IDX             = 0;
const unsigned STAT_INS_SIZE_MEAN_IDX      = 1;
const unsigned STAT_INS_SIZE_SD_IDX        = 2;
const unsigned STAT_INS_SIZE_MEDIAN_IDX    = 3;

const unsigned STAT_QUALITY_MEAN_IDX       = 4;
const unsigned STAT_QUALITY_SD_IDX         = 5;
const unsigned STAT_QUALITY_MEDIAN_IDX     = 6;

const unsigned STAT_SINGLES_MEAN_IDX       = 7;
const unsigned STAT_SINGLES_SD_IDX         = 8;
const unsigned STAT_SINGLES_MEDIAN_IDX     = 9;

const unsigned STAT_READ1_LEN_IDX          = 10;
const unsigned STAT_READ2_LEN_IDX          = 11;
const unsigned STAT_REL_ORIENT_IDX         = 12;


/* ----- ----- ----- ----- ----- -----
 *
 * ----- Auxiliary functions     -----
 *
 * ----- ----- ----- ----- ----- ----- */



bool
calcStats(std::vector<int32_t>& data,
          PairStatSet& stats) {

    const unsigned n(data.size());
    if (n == 0) {
        stats.clear();
        return false;
    }

    // First we'll sort to calculate the median
    sort(data.begin(), data.end());
    stats.median = data[(uint32_t)(n * 0.5)];

    // this is iqr, but we call it sd:
    stats.sd = data[(uint32_t)(n * 0.75)] - data[(uint32_t)(n * 0.25)];

    double sum(0);
    for (unsigned ii(0); ii < n; ii++) {
        sum += data[ii];
    }
    stats.mean = sum / (double)n;
    return true;
}



static
bool
isStatSetMatch(const PairStatSet& a,
               const PairStatSet& b) {

    static const double statsPrecision(0.005);

    if (std::abs(a.median-b.median)>=statsPrecision) return false;
    if (std::abs(a.sd-b.sd)>=statsPrecision) return false;
    if (std::abs(a.mean-b.mean)>=statsPrecision) return false;
    return true;
}



// This produces a useful result only when both reads align to the same
// chromosome.
static
ReadPairOrient
getRelOrient(const bam_record& br) {

    pos_t pos1 = br.pos();
    bool is_fwd_strand1 = br.is_fwd_strand();
    pos_t pos2 = br.mate_pos();
    bool is_fwd_strand2 = br.is_mate_fwd_strand();

    if (! br.is_first()) {
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
std::ostream&
operator<<(std::ostream& os, const PairStatSet& pss) {
    os << pss.mean << '\t' << pss.sd << '\t' << pss.median;
    return os;
}


/* ----- ----- ----- ----- ----- -----
 *
 * ----- ReadPairStats  -----
 *
 * ----- ----- ----- ----- ----- ----- */

ReadGroupStats::
ReadGroupStats(const std::vector<std::string>& data) {

    using namespace illumina::blt_util;

    // Initialize data
    InsSize.mean = parse_double_str(data[STAT_INS_SIZE_MEAN_IDX]);
    InsSize.sd = parse_double_str(data[STAT_INS_SIZE_SD_IDX]);
    InsSize.median = parse_double_str(data[STAT_INS_SIZE_MEDIAN_IDX]);

    // Limiting to 2 reads per fragment
    readLens.resize(2);
    readLens[0] = parse_int_str(data[STAT_READ1_LEN_IDX]);
    readLens[1] = parse_int_str(data[STAT_READ2_LEN_IDX]);

    relOrients.setVal(PAIR_ORIENT::get_index(data[STAT_REL_ORIENT_IDX].c_str()));
}



// set read pair statistics from a bam reader object:
//
ReadGroupStats::
ReadGroupStats(const std::string& bamFile) {

    static const unsigned statsCheckCnt(100000);
    static const unsigned maxPosCount(1);

    readLens.resize(2,0);

    bam_streamer read_stream(bamFile.c_str());

    const bam_header_t& header(* read_stream.get_header());
    const int32_t nChrom(header.n_targets);
    std::vector<int32_t> chromSize(nChrom,0);
    std::vector<int32_t> chromHighestPos(nChrom,-1);
    for (int32_t i(0); i<nChrom; ++i) {
        chromSize[i] = (header.target_len[i]);
    }

    bool isConverged(false);
    bool isStopEstimation(false);
    unsigned recordCnts(0);

    unsigned posCount(0);

    bool isPairTypeSet(false);
    bool isRead1Set(false);
    bool isRead2Set(false);

    PairStatsData psd;

    bool isActiveChrom(true);
    while (isActiveChrom && (!isStopEstimation)) {
        isActiveChrom=false;
        for (int32_t i(0); i<nChrom; ++i) {
            if (isStopEstimation) break;

            const int32_t startPos(chromHighestPos[i]+1);
#ifdef DEBUG_RPS
            std::cerr << "INFO: Stats requesting bam region starting from: chrid: " << i << " start: " << startPos << "\n";
#endif
            read_stream.set_new_region(i,startPos,chromSize[i]);
            while (read_stream.next()) {
                const bam_record& al(*(read_stream.get_record_ptr()));
                if (al.pos()<startPos) continue;

                if (al.pos()!=chromHighestPos[i]) {
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

                if (! isPairTypeSet) {
                    relOrients = getRelOrient(al);

                    if (readNum==1) isRead1Set=true;
                    if (readNum==2) isRead2Set=true;

                    readLens[readNum-1] = al.read_size();
                    isPairTypeSet=(isRead1Set && isRead2Set);
                }

                //  NOTE: Apparently we should not remove length from insert size...
                psd.fragmentLengths.push_back((std::abs(al.template_size()) - al.read_size()));

                if ((recordCnts % statsCheckCnt) != 0) continue;

#ifdef DEBUG_RPS
                log_os << "INFO: Checking stats convergence at record count : " << recordCnts << "'\n"
                       << "INFO: Stats before convergence check: ";
                store(log_os);
                log_os << "\n";
#endif

                isConverged=computePairStats(psd);
                if (isConverged || (recordCnts>5000000)) isStopEstimation=true;
                break;
            }
        }
    }

    if (! isConverged) {
        if (psd.fragmentLengths.size()<1000) {
            log_os << "ERROR: Can't generate pair statistics for BAM file " << bamFile << "\n";
            log_os << "\tTotal observed read pairs: " << psd.fragmentLengths.size() << "\n";
            exit(EXIT_FAILURE);
        }
        isConverged=computePairStats(psd);
    }
    if (! isConverged) {
        log_os << "WARNING: read pair statistics did not converge\n";

        // make sure stats are estimated for all values if we're going to continue:
        computePairStats(psd,true);
    }
}



unsigned
ReadGroupStats::
getReadLen(const unsigned readNum) const {
    assert(readNum>0);
    return readLens[readNum - 1];
}



bool
ReadGroupStats::
computePairStats(PairStatsData& psd, const bool isForcedConvergence) {

    // Calculate new mean, median and sd
    PairStatSet newVals;
    const bool calcStatus(calcStats(psd.fragmentLengths, newVals));
    if (! isForcedConvergence) {
        if (! calcStatus) return false;

        if (! isStatSetMatch(InsSize,newVals)) {
            InsSize = newVals;
            return false;
        }
    }

    return true;
}



void
ReadGroupStats::
store(std::ostream& os) const {
    os << InsSize << "\t"
       << readLens[0] << "\t"
       << readLens[1] << "\t"
       << relOrients;
}

