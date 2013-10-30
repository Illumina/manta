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

#pragma once

#include "blt_util/align_path.hh"
#include "blt_util/qscore.hh"

#include <boost/utility.hpp>


/// provide semi-aligned score metric for an aligned read
///
/// note: implemented as singleton class which pre-computes score components
struct ReadScorer : private boost::noncopyable
{
    /// return semi aligned read score given a seq-match/mismatch cigar and qualities
    ///
    /// not readLen is not strictly required but used as a safety check here
    static
    double
    getSemiAlignedMetric(
        const unsigned readLen,
        const ALIGNPATH::path_t& apath,
        const uint8_t* qual)
    {
        return get().getSemiAlignedMetricImpl(readLen, apath,qual);
    }

private:
    ReadScorer();

    /// return singleton instance of ReadScorer
    static const ReadScorer& get()
    {
        static const ReadScorer rs;
        return rs;
    }

    double
    getSemiAlignedMetricImpl(
        const unsigned readLen,
        const ALIGNPATH::path_t& apath,
        const uint8_t* qual) const;

    double
    getLogRatio(const int qual) const
    {
        qphred_cache::qscore_check(qual, "basecall quality");
        return _logpcorrectratio[qual];
    }

    enum { MAX_QSCORE  = qphred_cache::MAX_QSCORE };
    double _logpcorrectratio[MAX_QSCORE+1];
};
