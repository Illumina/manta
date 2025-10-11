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

#pragma once

#include "blt_util/qscore_cache.hpp"

/// Provides a variation on regular qscores by incorporating SNP probability
///
/// SNP probability limits the error if we have to explain the difference between a read and the reference,
/// however it should not be used in cases where the biological variation has already been accounted for in
/// the comparison.
///
struct qscore_snp {
  qscore_snp(const double snp_prob);

  double qphred_to_error_prob(const int qscore) const
  {
    qphred_cache::qscore_check(qscore, "basecall quality");
    return _q2p[qscore];
  }

  double qphred_to_ln_comp_error_prob(const int qscore) const
  {
    qphred_cache::qscore_check(qscore, "basecall quality");
    return _q2lncompe[qscore];
  }

  double qphred_to_ln_error_prob(const int qscore) const
  {
    qphred_cache::qscore_check(qscore, "basecall quality");
    return _q2lne[qscore];
  }

private:
  double _q2p[qphred_cache::MAX_QSCORE + 1];
  double _q2lncompe[qphred_cache::MAX_QSCORE + 1];
  double _q2lne[qphred_cache::MAX_QSCORE + 1];
};
