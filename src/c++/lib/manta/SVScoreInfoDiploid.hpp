//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Chris Saunders and Xiaoyu Chen
///

#pragma once

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <array>
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

namespace DIPLOID_GT {
enum index_t { REF, HET, HOM, SIZE };

/// Prior probability of expected alt allele for each genotype in the order of {REF, HET, HOM}
/// Note the alt prior for HOM is set to 0.99, allowing minor evidence (0.01) for ref allele
static const std::array<float, SIZE> altPriors = {0., 0.5, 0.99};
// pre-compute log values for cache
static const std::array<float, SIZE> altLnPriors = {
    std::log(altPriors[REF]), std::log(altPriors[HET]), std::log(altPriors[HOM])};
static const std::array<float, SIZE> altLnCompPriors = {
    std::log(1 - altPriors[REF]), std::log(1 - altPriors[HET]), std::log(1 - altPriors[HOM])};

inline const char* label(const index_t i)
{
  switch (i) {
  case REF:
    return "ref";
  case HET:
    return "het";
  case HOM:
    return "hom";
  default:
    assert(false && "Unknown GT state");
    return nullptr;
  }
}

inline const char* label(const unsigned i)
{
  return label(static_cast<index_t>(i));
}

inline float altFraction(const index_t i)
{
  assert((i < SIZE) && "Unknown GT state");
  return altPriors[i];
}

inline double altLnFraction(const index_t i)
{
  assert((i < SIZE) && "Unknown GT state");
  return altLnPriors[i];
}

inline double altLnCompFraction(const index_t i)
{
  assert((i < SIZE) && "Unknown GT state");
  return altLnCompPriors[i];
}

}  // namespace DIPLOID_GT

struct SVScoreInfoDiploidSample {
  SVScoreInfoDiploidSample() : phredLoghood(DIPLOID_GT::SIZE, 0), pprob(DIPLOID_GT::SIZE, 0) {}

  void clear()
  {
    filters.clear();
    gt      = DIPLOID_GT::REF;
    gtScore = 0;
    std::fill(phredLoghood.begin(), phredLoghood.end(), 0);
    std::fill(pprob.begin(), pprob.end(), 0);
  }

  std::set<std::string> filters;

  DIPLOID_GT::index_t gt = DIPLOID_GT::REF;

  unsigned gtScore = 0;  ///< quality score of genotype

  std::vector<unsigned> phredLoghood;
  std::vector<double>   pprob;
};

std::ostream& operator<<(std::ostream& os, const SVScoreInfoDiploidSample& sid);

/// consolidate all germline scoring results applied to an SV candidate
struct SVScoreInfoDiploid {
  void setSampleCount(const unsigned sampleCount) { samples.resize(sampleCount); }

  void clear()
  {
    filters.clear();
    altScore = 0;
    for (auto& sample : samples) {
      sample.clear();
    }
  }

  std::set<std::string> filters;

  /// Quality score indicating any non-reference state (regardless of specific genotype)
  unsigned altScore = 0;

  std::vector<SVScoreInfoDiploidSample> samples;
};

std::ostream& operator<<(std::ostream& os, const SVScoreInfoDiploid& sid);
