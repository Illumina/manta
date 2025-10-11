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
