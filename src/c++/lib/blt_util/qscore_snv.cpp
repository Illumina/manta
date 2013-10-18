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
/// \author Chris Saunders
///

#include "blt_util/qscore_snv.hh"
#include "blt_util/qscore.hh"
#include "blt_util/math_util.hh"


qscore_snv::
qscore_snv(
    const double snv_prob)
{
    static const int MAX_QSCORE(qphred_cache::MAX_QSCORE);

    const double comp_snv3(1.-(snv_prob/3.));

    for (int i(0); i<=MAX_QSCORE; ++i)
    {
        const double qerr(phred_to_error_prob(static_cast<double>(i)));
        _q2p[i] = (qerr * comp_snv3) + ((1-qerr) * snv_prob);
        _q2lncompe[i] = log1p_switch(-_q2p[i]);
        _q2lne[i] = std::log(_q2p[i]);
    }
}
