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

#pragma once

#include <string>
#include <vector>

#include "EdgeOptions.hpp"
#include "common/Program.hpp"
#include "options/AlignmentFileOptions.hpp"
#include "options/CallOptionsDiploid.hpp"
#include "options/CallOptionsShared.hpp"
#include "options/CallOptionsSomatic.hpp"
#include "options/CallOptionsTumor.hpp"
#include "options/ReadScannerOptions.hpp"
#include "options/SVRefinerOptions.hpp"

struct GSCOptions {
  bool isSomatic() const { return (!somaticOutputFilename.empty()); }

  bool isTumorOnly() const { return (!tumorOutputFilename.empty()); }

  bool isGenerateEvidenceBam() const { return (!evidenceBamStub.empty()); }

  AlignmentFileOptions alignFileOpt;
  EdgeOptions          edgeOpt;
  ReadScannerOptions   scanOpt;
  SVRefinerOptions     refineOpt;
  CallOptionsShared    callOpt;
  CallOptionsDiploid   diploidOpt;
  CallOptionsSomatic   somaticOpt;
  CallOptionsTumor     tumorOpt;

  int workerThreadCount = 1;

  std::string graphFilename;
  std::string referenceFilename;
  std::string statsFilename;
  std::string chromDepthFilename;
  std::string edgeRuntimeFilename;
  std::string edgeStatsFilename;
  std::string edgeStatsReportFilename;

  std::string candidateOutputFilename;
  std::string diploidOutputFilename;
  std::string somaticOutputFilename;
  std::string tumorOutputFilename;
  std::string rnaOutputFilename;
  std::string evidenceBamStub;

  /// Provide some high-level log info to assist in debugging
  bool isVerbose = false;

  /// If true, skip assembly and run a low-resolution, breakdancer-like subset of the workflow
  bool isSkipAssembly = false;

  /// If true, skip quality scoring and output candidates only
  bool isSkipScoring = false;

  /// If true, turn on retrieval of poorly mapped remote reads for assembly
  bool enableRemoteReadRetrieval = false;

  /// If true, RNA specific filtering on candidates and diploid scoring is used
  bool isRNA = false;

  /// For unstranded RNA data, the direction of fusion transcripts is unknown
  bool isUnstrandedRNA = false;

  /// How many spanning evidence observations are required to become a candidate?
  unsigned minCandidateSpanningCount = 3;

  /// Min size for scoring and scored output following candidate generation
  unsigned minScoredVariantSize = 50;

  /// if true, an assembled contig is written in VCF
  bool isOutputContig = false;

  /// if true, turn off the filter of insignificant evidence signal
  bool skipEvidenceSignalFilter = false;
};

void parseGSCOptions(const illumina::Program& prog, int argc, char* argv[], GSCOptions& opt);
