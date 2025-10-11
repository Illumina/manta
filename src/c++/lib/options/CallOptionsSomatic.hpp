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

#pragma once

#include "boost/program_options.hpp"

struct CallOptionsSomatic {
  float germlineSVPrior = 1e-5f;
  float somaticSVPrior  = 1e-7f;

  /// small/large values below reflect our expectation that there is more shared small event noise in small
  /// events
  ///
  /// expected shared tumor-normal sample noise rates for "small" SVs, ramp is from 3k->5k for small to large.
  float smallNoiseSVPrior = 1e-9f;
  /// expected shared tumor-normal sample noise rates for "large" SVs
  float largeNoiseSVPrior = 1e-10f;

  /// breakpoints where the non-tumor depth is greater than the chromosome average x this factor
  /// are filtered out:
  float       maxDepthFactor      = 3.0f;
  std::string maxDepthFilterLabel = "MaxDepth";

  /// minimum somatic quality to print out a somatic variant
  unsigned minOutputSomaticScore = 10;
  /// minimum somatic quality which passes vcf filtration
  unsigned    minPassSomaticScore  = 30;
  std::string minSomaticScoreLabel = "MinSomaticScore";

  float       maxMQ0Frac      = 0.4f;
  std::string maxMQ0FracLabel = "MaxMQ0Frac";
};

boost::program_options::options_description getOptionsDescription(CallOptionsSomatic& opt);
