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
/// \author Chris Saunders
///

#include "manta/SVScoreInfoRna.hpp"
#include "blt_util/log.hpp"

#include <iostream>

const std::string SVScoreInfoRna::rnaFilterLabel = "LowEvidence";
const std::string SVScoreInfoRna::impreciseLabel = "Imprecise";
const std::string SVScoreInfoRna::localLabel     = "Local";

std::ostream& operator<<(std::ostream& os, const SVScoreInfoRna& sid)
{
  os << "RnaSVScoreInfo "
     << " altScore=" << sid.altScore;
  os << " filters:";
  for (const std::string& filter : sid.filters) {
    os << " " << filter;
  }
  os << "\n";
  return os;
}
