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

#include "ESLOptions.hpp"
#include "svgraph/SVLocusSet.hpp"

#include <memory>
#include <string>

struct bam_streamer;

/// Provides SV loci estimation methods over multiple regions, and manages the one-time initialization costs
/// of this process
struct EstimateSVLociRunner : private boost::noncopyable {
  /// \param[in] opt Options for estimation process
  explicit EstimateSVLociRunner(const ESLOptions& opt);

  /// Run the SVlocus estimation process and the specified region and merge results into \p mergedSet
  ///
  /// \param[in] region Target region for estimation process
  void estimateSVLociForSingleRegion(const std::string& region);

  /// \brief Provide const access to the SV locus graph that this object is building.
  const SVLocusSet& getLocusSet() const { return *_mergedSetPtr; }

private:
  const ESLOptions _opt;

  std::vector<std::shared_ptr<bam_streamer>> _bamStreams;

  /// Estimated SVlocus graph components should be merged into this object
  std::shared_ptr<SVLocusSet> _mergedSetPtr;
};
