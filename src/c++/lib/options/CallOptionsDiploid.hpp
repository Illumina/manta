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

#pragma once

#include "boost/program_options.hpp"

struct CallOptionsDiploid {
  float indelPrior = 1e-5f;

  /// breakpoints where the non-tumor depth is greater than the chromosome average x this factor are filtered
  /// out:
  float       maxDepthFactor      = 3.0f;
  std::string maxDepthFilterLabel = "MaxDepth";

  /// minimum QUAL score to print out a diploid variant
  unsigned minOutputAltScore = 10;

  /// minimum QUAL score to PASS a diploid variant
  unsigned    minPassAltScore   = 20;
  std::string minAltFilterLabel = "MinQUAL";

  // control filtration based on MQ0 fraction:
  float       maxMQ0Frac      = 0.4f;
  std::string maxMQ0FracLabel = "MaxMQ0Frac";

  /// filter for large SVs with no pair support
  std::string noPairSupportLabel = "NoPairSupport";

  /// no sample passes all sample-specific filters
  std::string failedSampleFTLabel = "SampleFT";

  /// below this GQ value, the SAMPLE filter is marked in the VCF
  unsigned    minPassGTScore   = 15;
  std::string minGTFilterLabel = "MinGQ";

  /// the SAMPLE filter for home-ref calls
  std::string homRefLabel = "HomRef";
};

boost::program_options::options_description getOptionsDescription(CallOptionsDiploid& opt);
