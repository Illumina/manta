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

#include "blt_util/qscore_cache.hpp"
#include "blt_util/blt_exception.hpp"
#include "blt_util/math_util.hpp"
#include "blt_util/qscore.hpp"

#include <iostream>
#include <sstream>

qphred_cache::qphred_cache()
{
  static const double q2lnp(-std::log(10.) / 10.);

  for (int i(0); i <= MAX_QSCORE; ++i) {
    q2p[i]       = phred_to_error_prob(static_cast<double>(i));
    q2lncompe[i] = log1p_switch(-q2p[i]);
    q2lne[i]     = static_cast<double>(i) * q2lnp;
    for (int j(0); j <= MAX_MAP; ++j) {
      mappedq[j][i] = error_prob_to_qphred(phred_to_mapped_error_prob(i, j));
    }
  }
}

void qphred_cache::invalid_qscore_error(const int qscore, const char* label)
{
  std::stringstream oss;
  oss << "Attempting to lookup invalid " << label << " score: " << qscore;
  throw blt_exception(oss.str().c_str());
}

void qphred_cache::high_qscore_error(const int qscore, const char* label)
{
  std::stringstream oss;
  oss << "Attempting to lookup " << label << " score " << qscore << " which exceeds the maximum cached "
      << label << " score of " << MAX_QSCORE;
  throw blt_exception(oss.str().c_str());
}
