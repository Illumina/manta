// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/qscore_cache.hh"


/// this object provides a variation on regular qscores by incorporating SNP probability
///
/// SNV probability limits the error if we have to explain the difference
/// between a read and the reference, however it should not be used in cases where the
/// biological variation has already been accounted for in the comparison.
///
struct qscore_snp
{
    qscore_snp(const double snp_prob);

    double
    qphred_to_error_prob(const int qscore) const
    {
        qphred_cache::qscore_check(qscore, "basecall quality");
        return _q2p[qscore];
    }

    double
    qphred_to_ln_comp_error_prob(const int qscore) const
    {
        qphred_cache::qscore_check(qscore, "basecall quality");
        return _q2lncompe[qscore];
    }

    double
    qphred_to_ln_error_prob(const int qscore) const
    {
        qphred_cache::qscore_check(qscore, "basecall quality");
        return _q2lne[qscore];
    }

private:
    double _q2p[qphred_cache::MAX_QSCORE+1];
    double _q2lncompe[qphred_cache::MAX_QSCORE+1];
    double _q2lne[qphred_cache::MAX_QSCORE+1];
};
