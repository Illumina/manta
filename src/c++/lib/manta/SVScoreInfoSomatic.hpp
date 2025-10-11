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
/// \author Chris Saunders and Xiaoyu Chen
///

#pragma once

#include <cassert>
#include <cstdlib>

#include <cmath>
#include <iosfwd>
#include <set>
#include <string>

namespace SOMATIC_GT {
// TODO: estimated from tumor
const double SOMATIC_MUTATION_FREQ = 0.6;

enum index_t { REF, HET, HOM, SOM, NOISE, SIZE };

inline const char* label(const index_t i)
{
  switch (i) {
  case REF:
    return "ref";
  case HET:
    return "het";
  case HOM:
    return "hom";
  case SOM:
    return "som";
  case NOISE:
    return "noise";
  default:
    assert(false && "Unknown GT state");
    return nullptr;
  }
}

inline const char* label(const unsigned i)
{
  return label(static_cast<index_t>(i));
}

inline float altFraction(const index_t i, const float somaticFreq, const float noiseFreq)
{
  switch (i) {
  case REF:
    return 0;
  case HET:
    return 0.5;
  case HOM:
    return 1.0;
  case SOM:
    return somaticFreq;
  case NOISE:
    return noiseFreq;
  default:
    assert(false && "Unknown GT state");
    return 0;
  }
}

inline double altLnFraction(const index_t i, const double somaticFreq, const double noiseFreq)
{
  static const double val[] = {std::log(0.), std::log(0.5), std::log(1.)};

  switch (i) {
  case REF:
    return val[0];
  case HET:
    return val[1];
  case HOM:
    return val[2];
  case SOM:
    return std::log(somaticFreq);
  case NOISE:
    return std::log(noiseFreq);
  default:
    assert(false && "Unknown GT state");
    return 0;
  }
}

inline double altLnCompFraction(const index_t i, const double somaticFreq, const double noiseFreq)
{
  static const double val[] = {std::log(1.), std::log(0.5), std::log(0.)};

  switch (i) {
  case REF:
    return val[0];
  case HET:
    return val[1];
  case HOM:
    return val[2];
  case SOM:
    return std::log(1 - somaticFreq);
  case NOISE:
    return std::log(1 - noiseFreq);
  default:
    assert(false && "Unknown GT state");
    return 0;
  }
}
}  // namespace SOMATIC_GT

/// consolidate all somatic scoring results applied to an SV candidate
struct SVScoreInfoSomatic {
  void clear()
  {
    filters.clear();
    somaticScore     = 0;
    somaticScoreTier = 0;
  }

  std::set<std::string> filters;

  unsigned      somaticScore     = 0;
  unsigned char somaticScoreTier = 0;
};

std::ostream& operator<<(std::ostream& os, const SVScoreInfoSomatic& sis);
