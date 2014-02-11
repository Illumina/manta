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
/// \author Chris Saunders
///

#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "blt_util/math_util.hh"
#include "blt_util/qscore.hh"
#include "blt_util/qscore_cache.hh"

#include <cstdlib>
#include <iostream>
#include <sstream>



qphred_cache::
qphred_cache()
{
    static const double q2lnp(-std::log(10.)/10.);

    for (int i(0); i<=MAX_QSCORE; ++i)
    {
        q2p[i] = phred_to_error_prob(static_cast<double>(i));
        q2lncompe[i] = log1p_switch(-q2p[i]);
        q2lne[i] = static_cast<double>(i)*q2lnp;
        for (int j(0); j<=MAX_MAP; ++j)
        {
            mappedq[j][i] = error_prob_to_qphred(phred_to_mapped_error_prob(i,j));
        }
    }
}



void
qphred_cache::
invalid_qscore_error(const int qscore,
                     const char* label)
{
    std::stringstream oss;
    oss << "ERROR: Attempting to lookup invalid " << label << " score: " << qscore;
    throw blt_exception(oss.str().c_str());
}


void
qphred_cache::
high_qscore_error(const int qscore,
                  const char* label)
{
    std::stringstream oss;
    oss << "ERROR: Attempting to lookup " << label << " score " << qscore << " which exceeds the maximum cached " << label << " score of " <<  MAX_QSCORE;
    throw blt_exception(oss.str().c_str());
}



qlogodds_cache::
qlogodds_cache() : q2p(q2p_base-MIN_QSCORE)
{
    for (int i(MIN_QSCORE); i<=MAX_QSCORE; ++i)
    {
        q2p[i] = logodds_to_error_prob(static_cast<double>(i));
    }
}



void
qlogodds_cache::
qscore_error(const int qscore) const
{
    log_os << "ERROR:: invalid logodds qscore: " << qscore << "\n";
    exit(EXIT_FAILURE);
}

