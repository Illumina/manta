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

#include <string>
#include <vector>

/// \brief Input alignment file object shared by all programs which require these as input
struct AlignmentFileOptions {
  unsigned diploidSampleCount() const
  {
    unsigned count(0);
    for (const bool isTumor : isAlignmentTumor) {
      if (!isTumor) count++;
    }
    return count;
  }

  typedef std::vector<std::string> files_t;

  files_t alignmentFilenames;

  /// Indicates which positions in the alignmnetFilename correspond to tumor
  std::vector<bool> isAlignmentTumor;
};
