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

struct CallOptionsTumor {
  /// breakpoints where the non-tumor depth is greater than the chromosome average x this factor
  /// are filtered out:
  float       maxDepthFactor      = 3.0f;
  std::string maxDepthFilterLabel = "MaxDepth";

  float       maxMQ0Frac      = 0.4f;
  std::string maxMQ0FracLabel = "MaxMQ0Frac";
};

boost::program_options::options_description getOptionsDescription(CallOptionsTumor& opt);
