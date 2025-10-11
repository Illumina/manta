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

#include <cassert>
#include <cstdint>

#include <array>

/// \brief Helper class used to precompute/cache quality scores
///
/// This is a singleton helper class used to accelerate quality services provided by qscore.hh.
struct qphred_cache {
  static double get_error_prob(const int qscore) { return qc().get_error_prob_imp(qscore); }

  static double get_ln_comp_error_prob(const int qscore) { return qc().get_lncompe_imp(qscore); }

  static double get_ln_error_prob(const int qscore) { return qc().get_lne_imp(qscore); }

  static int get_mapped_qscore(const int basecall_qscore, const int mapping_qscore)
  {
    return qc().get_mapped_qscore_imp(basecall_qscore, mapping_qscore);
  }

  enum { MAX_QSCORE = 70, MAX_MAP = 90 };

  static void qscore_check(const int qscore, const char* label)
  {
    if (qscore < 0) invalid_qscore_error(qscore, label);
    if (qscore > MAX_QSCORE) high_qscore_error(qscore, label);
  }

private:
  qphred_cache();

  static const qphred_cache& qc()
  {
    static const qphred_cache qc;
    return qc;
  }

  static void invalid_qscore_error(const int qscore, const char* label);
  static void high_qscore_error(const int qscore, const char* label);

  static void qscore_check_int(const int qscore)
  {
    static const char* label = "phred";
    qscore_check(qscore, label);
  }

  double get_error_prob_imp(const int qscore) const
  {
    qscore_check_int(qscore);
    return q2p[qscore];
  }

  double get_lncompe_imp(const int qscore) const
  {
    qscore_check_int(qscore);
    return q2lncompe[qscore];
  }

  double get_lne_imp(const int qscore) const
  {
    qscore_check_int(qscore);
    return q2lne[qscore];
  }

  int get_mapped_qscore_imp(const int basecall_qscore, int mapping_qscore) const
  {
    static const char* label = "basecall quality";
    qscore_check(basecall_qscore, label);
    assert(mapping_qscore >= 0);
    if (mapping_qscore > MAX_MAP) {
      mapping_qscore = MAX_MAP;
    }
    return mappedq[mapping_qscore][basecall_qscore];
  }

  std::array<double, MAX_QSCORE + 1> q2p;
  std::array<double, MAX_QSCORE + 1> q2lncompe;
  std::array<double, MAX_QSCORE + 1> q2lne;
  uint8_t                            mappedq[MAX_MAP + 1][MAX_QSCORE + 1];
};
