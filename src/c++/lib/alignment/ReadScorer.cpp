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
/// \author Ole Schulz-Trieglaff
///

#include "ReadScorer.hh"
#include "blt_util/math_util.hh"

#include "boost/foreach.hpp"

#include <cassert>

//#define DEBUG_RS


#ifdef DEBUG_RS
#include "blt_util/log.hh"
#include <iostream>
#endif



ReadScorer::
ReadScorer()
{
#ifdef DEBUG_RS
    log_os << "Filling logpcorrectratio table" << std::endl;
#endif
    // below we take log(1+(-i)) so set a dummy value for i = 0
    assert(MAX_QSCORE>=0);
    _logpcorrectratio[0] = 0.99;
    for (unsigned i(1); i<=MAX_QSCORE; ++i)
    {
        const double eprob(qphred_to_error_prob(i));
#ifdef DEBUG_RS
        log_os << "i=" << i << " " << eprob << " " << (log1p_switch(-eprob) - std::log(eprob/3.)) << std::endl;
#endif
        _logpcorrectratio[i] = log1p_switch(-eprob) - std::log(eprob/3);
    }
}



double
ReadScorer::
getSemiAlignedMetricImpl(
    const unsigned readLen,
    const ALIGNPATH::path_t& apath,
    const uint8_t* qual) const
{
    using namespace ALIGNPATH;

    unsigned posInRead = 0;
    double score(0.);

    BOOST_FOREACH(const path_segment& ps, apath)
    {
        assert((ps.type != MATCH) && "Incorrect CIGAR type, matches must be converted to SEQ_MATCH/SEQ_MISMATCH");

        if ((ps.type==SOFT_CLIP) || (ps.type==SEQ_MISMATCH))
        {
            assert((posInRead+ps.length) <= readLen);
            for (unsigned j(0); j<ps.length; ++j)
            {
                score += getLogRatio(qual[posInRead+j]);
            }
        }
        if (is_segment_type_read_length(ps.type)) posInRead += ps.length;
    }

    return score;
}
