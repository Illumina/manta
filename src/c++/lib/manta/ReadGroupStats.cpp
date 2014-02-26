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
#include "common/Exceptions.hh"
#include "manta/SVLocusScanner.hh"

#include "boost/foreach.hpp"

#include <iostream>
#include <set>
#include <sstream>
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
PAIR_ORIENT::index_t
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

    return PAIR_ORIENT::get_index(pos1,is_fwd_strand1,pos2,is_fwd_strand2);
}



struct ReadGroupTmpInfo
{
    explicit
    ReadGroupTmpInfo(const char* rgLabel = NULL) :
        _isFinalized(false),
        _totalOrientCount(0),
        _rgLabel(NULL==rgLabel ? "" : rgLabel)
    {
        std::fill(_orientCount.begin(),_orientCount.end(),0);
    }

    void
    addOrient(
        const bam_record& bamRead)
    {
        if (bamRead.pos() == bamRead.mate_pos()) return;

        static const unsigned maxOrientCount(100000);
        if (getOrientCount() < maxOrientCount)
        {
            addOrient(getRelOrient(bamRead));
        }
    }

    void
    addOrient(
        const PAIR_ORIENT::index_t id)
    {
        assert(! _isFinalized);
        assert((id>=0) && (id<PAIR_ORIENT::SIZE));

        _orientCount[id]++;
        _totalOrientCount++;
    }

    unsigned
    getOrientCount() const
    {
        return _totalOrientCount;
    }

    const ReadPairOrient&
    getOrient()
    {
        finalize();
        return _finalOrient;
    }

private:
    void
    finalize()
    {
        if (_isFinalized) return;

        bool isMaxIndex(false);
        unsigned maxIndex(0);
        for(unsigned i(0);i<_orientCount.size();++i)
        {
            if ((! isMaxIndex) || (_orientCount[i] > _orientCount[maxIndex]))
            {
                isMaxIndex=true;
                maxIndex=i;
            }
        }

        assert(isMaxIndex);

        _finalOrient.setVal(maxIndex);

        { // make sure there's a dominant consensus orientation and that we have a minimum number of samples:
            static const unsigned minCount(100);
            static const float minMaxFrac(0.9);

            using namespace illumina::common;

            if (_totalOrientCount < minCount)
            {
                std::ostringstream oss;
                oss << "ERROR: Too few observations (" << _totalOrientCount << ") to determine pair orientation for read group '" << _rgLabel << "'\n";
                BOOST_THROW_EXCEPTION(LogicException(oss.str()));
            }

            const unsigned minMaxCount(minMaxFrac*_totalOrientCount);
            if (_orientCount[maxIndex] < minMaxCount)
            {
                const unsigned maxPercent((_orientCount[maxIndex]*100)/_totalOrientCount);
                std::ostringstream oss;
                oss << "ERROR: Can't determine consensus pair orientation of read group '" << _rgLabel << "'.\n"
                    << "\tMost frequent orientation is '" << _finalOrient << "' (" << maxPercent << "% of " << _totalOrientCount << " total observations)\n";
                BOOST_THROW_EXCEPTION(LogicException(oss.str()));
            }
        }

        _isFinalized=true;
    }

    bool _isFinalized;
    unsigned _totalOrientCount;
    std::string _rgLabel;
    boost::array<unsigned,PAIR_ORIENT::SIZE> _orientCount;

    ReadPairOrient _finalOrient;
};



struct ReadAlignFilter
{
    /// use only the most conservative alignments to generate fragment stats --
    /// filter reads containing any cigar types besides MATCH:
    bool
    isFilterRead(
        const bam_record& bamRead)
    {
        using namespace ALIGNPATH;

        bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), _apath);

        BOOST_FOREACH(const path_segment& ps, _apath)
        {
            if (! is_segment_align_match(ps.type)) return true;
        }
        return false;
    }

private:
    ALIGNPATH::path_t _apath;
};



/// this object handles filtration of reads which are:
///
/// 1. not downstream or (if both reads start at same position) not the second in order in the bam file
/// 2. part of a high depth pileup region
///
struct ReadPairDepthFilter
{
    ReadPairDepthFilter() :
        _posCount(0),
        _lastTargetId(0),
        _lastPos(0)
    {}

    bool
    isFilterRead(
        const bam_record& bamRead)
    {
        static const unsigned maxPosCount(1);

        if (bamRead.target_id() != _lastTargetId)
        {
            _goodMates.clear();
            _lastTargetId=bamRead.target_id();
            _posCount=0;
            _lastPos = bamRead.pos();
        }
        else if(bamRead.pos() != _lastPos)
        {
            _posCount=0;
            _lastPos = bamRead.pos();
        }

        // sample each read pair once by sampling stats from
        // downstream read only, or whichever read is encountered
        // second if the read and its mate start at the same position:
        const bool isDownstream(bamRead.pos() > bamRead.mate_pos());
        const bool isSamePos(bamRead.pos() == bamRead.mate_pos());

        if (isDownstream || isSamePos)
        {
            const int mateReadNo( bamRead.is_first() ? 2 : 1);
            const ReadKey mateKey(bamRead.qname(), mateReadNo);

            mateMap_t::iterator i(_goodMates.find(mateKey));

            if (i == _goodMates.end())
            {
                if (isDownstream) return true;
            }
            else
            {
                _goodMates.erase(i);
                return false;
            }
        }

        // to prevent high-depth pileups from overly biasing the
        // read stats, we only take maxPosCount read pairs from each start
        // pos. by not inserting a key in goodMates, we also filter
        // the downstream mate:
        if (_posCount>=maxPosCount) return true;
        ++_posCount;

        /// crude mechanism to manage total set memory
        static const unsigned maxMateSetSize(100000);
        if (_goodMates.size() > maxMateSetSize) _goodMates.clear();

        _goodMates.insert(ReadKey(bamRead));

        return true;
    }

private:
    typedef std::set<ReadKey> mateMap_t;

    unsigned _posCount;
    mateMap_t _goodMates;

    int _lastTargetId;
    int _lastPos;
};



// set read pair statistics from a bam reader object:
//
ReadGroupStats::
ReadGroupStats(
    const std::string& statsBamFile)
{
    static const unsigned statsCheckCnt(100000);

    bam_streamer read_stream(statsBamFile.c_str());

    const bam_header_t& header(* read_stream.get_header());
    const int32_t chromCount(header.n_targets);
    std::vector<int32_t> chromSize(chromCount,0);
    std::vector<int32_t> chromHighestPos(chromCount,-1);
    for (int32_t i(0); i<chromCount; ++i)
    {
        chromSize[i] = (header.target_len[i]);
    }

    bool isConverged(false);
    bool isStopEstimation(false);
    bool isFirstEstimation(true);
    SizeDistribution oldFragSize;

    ReadAlignFilter alignFilter;
    ReadPairDepthFilter pairFilter;
    unsigned recordCount(0);
    bool isActiveChrom(true);

    ReadGroupTmpInfo rgInfo(statsBamFile.c_str());

    while (isActiveChrom && (!isStopEstimation))
    {
        isActiveChrom=false;
        for (int32_t chromIndex(0); chromIndex<chromCount; ++chromIndex)
        {
            if (isStopEstimation) break;

            const int32_t startPos(chromHighestPos[chromIndex]+1);
#ifdef DEBUG_RPS
            std::cerr << "INFO: Stats requesting bam region starting from: chrid: " << chromIndex << " start: " << startPos << "\n";
#endif
            read_stream.set_new_region(chromIndex,startPos,chromSize[chromIndex]);
            while (read_stream.next())
            {
                const bam_record& bamRead(*(read_stream.get_record_ptr()));
                if (bamRead.pos()<startPos) continue;

                chromHighestPos[chromIndex]=bamRead.pos();
                isActiveChrom=true;

                // filter common categories of undesirable reads:
                if (SVLocusScanner::isReadFilteredCore(bamRead)) continue;

                if (! is_mapped_chrom_pair(bamRead)) continue;
                if (bamRead.map_qual()==0) continue;

                // filter any split reads with an SA tag:
                static const char SAtag[] = {'S','A'};
                if (NULL != bamRead.get_string_tag(SAtag)) continue;

                // remove reads without perfect alignments
                if (alignFilter.isFilterRead(bamRead)) continue;

                /// filter out upstream reads and high depth regions:
                if (pairFilter.isFilterRead(bamRead)) continue;

                // Assert only two reads per fragment
                const unsigned readNum(bamRead.is_first() ? 1 : 2);
                assert(bamRead.is_second() == (readNum == 2));

                // get orientation stats before final filter for innie reads below:
                //
                // we won't use anything but innie reads, but sampling up hear allows us to detect, ie. a mate-pair library
                // so that we can blow-up with an informative error
                //
                rgInfo.addOrient(bamRead);

                // filter mapped innies on the same chrom
                //
                // note we don't rely on the proper pair bit because this already contains an arbitrary length filter
                // AND subjects the method to aligner specific variation
                //
                // TODO: ..note this locks-in standard ilmn orientation -- okay for now but this function needs major
                // re-arrangement for mate-pair support, we could still keep independence from each aligner's proper
                // pair decisions by estimating a fragment distro for each orientation and only keeping the one with
                // the most samples
                //
                if (! is_innie_pair(bamRead)) continue;

                // made it through all filters!
                ++recordCount;


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
            using namespace illumina::common;

            std::ostringstream oss;
            oss << "ERROR: Can't generate pair statistics for BAM file " << statsBamFile << "\n"
                << "\tTotal observed read pairs: " << fragStats.totalObservations() << "\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
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

    // finalize pair orientation:
    relOrients = rgInfo.getOrient();

    if (relOrients.val() != PAIR_ORIENT::Rp)
    {
        using namespace illumina::common;

        std::ostringstream oss;
        oss << "ERROR: Unexpected consensus read orientation (" << relOrients << ") for read group: " << statsBamFile << "\n"
            << "\tManta currently handles paired-end (FR) reads only.\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    // final step before saving is to cut-off the extreme end of the fragment size distribution, this
    // is similar the some aligner's proper-pair bit definition of (3x the standard mean, etc.)
    static const float filterQuant(0.9995);
    fragStats.filterObservationsOverQuantile(filterQuant);
}
