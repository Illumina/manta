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


#include "SizeDistribution.hh"

#include "blt_util/log.hh"

#include "boost/foreach.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

//#define DEBUG_RPS

#if 0
static
void
writeFragSizeHashItem(std::ostream& os, const PairStatSet::frag_map_type& _fragMap, const int k)
{
    os << k << " -> <"
       << _fragMap.find(k)->second.first
       << ", "
       << _fragMap.find(k)->second.second
       << ">\n";
}
#endif



static
void
populateCdfQuantiles(
    SizeDistribution::map_type& sizeMap,
    const unsigned totalCount,
    std::vector<float>& quantiles)
{
    const unsigned quantileNum(quantiles.size());
    const float pFactor(1./static_cast<float>(totalCount));

    unsigned fillBase(0);
    unsigned cumulativeCount(0);
    BOOST_REVERSE_FOREACH(SizeDistribution::map_type::value_type& val, sizeMap)
    {
        cumulativeCount += (val.second.count);
        assert(cumulativeCount <= totalCount);

        // update the hash map with cumulative prob value
        val.second.cprob = (cumulativeCount * pFactor);

        const unsigned fillNext = rint(val.second.cprob * quantileNum);
        for (; fillBase < fillNext; fillBase++)
        {
            quantiles[fillBase] = val.first;
        }
    }
}




void
SizeDistribution::
calcStats() const
{
#ifdef DEBUG_RPS
    log_os << "Calculating stats...\n"
           << "numOfSized=" << _sizeMap.size() << "\n";
#endif
    _isStatsComputed=true;
    if (_sizeMap.empty()) return;

    populateCdfQuantiles(_sizeMap, _totalCount, _quantiles);
}



int
SizeDistribution::
quantile(const float p) const
{
    assert((p >= 0.) && (p <= 1.));

    static const int maxBin(_quantileNum - 1);
    if(! _isStatsComputed) calcStats();

    int bin(static_cast<int>(ceil(p * _quantileNum) - 1));
    if(bin < 0) bin=0;
    if(bin > maxBin) bin=maxBin;
    return _quantiles[bin];
}



float
SizeDistribution::
cdf(const int size) const
{
    if(! _isStatsComputed) calcStats();

    const map_type::const_iterator sizeIter(_sizeMap.lower_bound(size));
    if(sizeIter == _sizeMap.end()) return 0;
    return sizeIter->second.cprob;
}



std::ostream&
operator<<(std::ostream& os, const SizeDistribution& sd)
{
    os << sd.totalObservations() << '\n';
    return os;
}
