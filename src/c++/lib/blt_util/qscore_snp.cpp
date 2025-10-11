//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
//
//

/// \file
/// \author Chris Saunders
///

#include "blt_util/qscore_snp.hpp"

#include "blt_util/math_util.hpp"
#include "blt_util/qscore.hpp"

qscore_snp::qscore_snp(const double snp_prob)
{
  static const int MAX_QSCORE(qphred_cache::MAX_QSCORE);

  const double comp_snp3(1. - (snp_prob / 3.));

  for (int i(0); i <= MAX_QSCORE; ++i) {
    const double qerr(phred_to_error_prob(static_cast<double>(i)));
    _q2p[i]       = (qerr * comp_snp3) + ((1 - qerr) * snp_prob);
    _q2lncompe[i] = log1p_switch(-_q2p[i]);
    _q2lne[i]     = std::log(_q2p[i]);
  }
}
