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

#include "alignment/AlignmentScores.hpp"
#include "options/IterativeAssemblerOptions.hpp"

typedef IterativeAssemblerOptions AssemblerOptions;

namespace ALIGNERPARAM {
static const int largeGapOpenScore(-100);
}

/// Options for the SV refiner step
///
/// Note that we have two categories of options for assembly and alignment,
/// one for small events, and one for large events
///
struct SVRefinerOptions {
  /// match, mismatch, open score ratios taken from bwa defaults (but not extend!) :
  ///
  SVRefinerOptions()
    : largeSVAlignScores(2, -8, -24, -1, -1),
      largeInsertEdgeAlignScores(2, -8, -18, -1, -1),
      largeInsertCompleteAlignScores(2, -8, ALIGNERPARAM::largeGapOpenScore, 0, -1),
      spanningAlignScores(2, -8, -12, -1, -1),
      largeGapOpenScore(ALIGNERPARAM::largeGapOpenScore),
      jumpScore(-100),
      RNAspanningAlignScores(2, -8, -19, -1, -1),
      RNAJumpScore(-100),
      RNAIntronOpenScore(-15),
      RNAIntronOffEdgeScore(-1),
      contigFilterScores(2, -8, -18, 0, -1)
  {
    // For breakend-spanning assemblies we require a larger contig than for small-variant assemblies
    spanningAssembleOpt.minContigLength = 75;

    // For breakend-spanning assemblies we require a larger contig than for small-variant assemblies
    RNAspanningAssembleOpt.minContigLength = 75;

    // Use smaller kmer for RNA
    RNAspanningAssembleOpt.minWordLength = 31;
  }

  // parameters for small SV assembly/alignment:
  AlignmentScores<int> largeSVAlignScores;  // large SV but at a single assembly locus
  AlignmentScores<int> largeInsertEdgeAlignScores;
  AlignmentScores<int> largeInsertCompleteAlignScores;
  AssemblerOptions     smallSVAssembleOpt;

  // parameters for large SV assembly/alignment:
  AlignmentScores<int> spanningAlignScores;
  const int            largeGapOpenScore;
  const int            jumpScore;
  AlignmentScores<int> RNAspanningAlignScores;
  const int            RNAJumpScore;
  const int            RNAIntronOpenScore;
  const int            RNAIntronOffEdgeScore;
  AssemblerOptions     spanningAssembleOpt;
  AssemblerOptions     RNAspanningAssembleOpt;

  /// Scores used to check for low quality assembly contig alignments
  AlignmentScores<int> contigFilterScores;
};
