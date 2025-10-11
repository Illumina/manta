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

///
/// \file
/// \author Trevor Ramsay
///

#include "testSVLocusScanner.hpp"

#include "test/testFileMakers.hpp"

#include "boost/make_unique.hpp"

std::unique_ptr<SVLocusScanner> buildTestSVLocusScanner(
    const bam_header_info& bamHeaderInfo, bool isRNA, int minCandidateVariantSizeInput)
{
  ReadScannerOptions opts      = ReadScannerOptions();
  opts.minCandidateVariantSize = minCandidateVariantSizeInput;
  opts.useOverlapPairEvidence  = isRNA;
  const ReadScannerOptions& constRefOpts(opts);

  TestStatsFileMaker             statsFile;
  TestAlignHeaderFileMaker       alignFile(bamHeaderInfo);
  const std::vector<std::string> alignFilenameVector = {alignFile.getFilename()};

  return boost::make_unique<SVLocusScanner>(
      constRefOpts, statsFile.getFilename(), alignFilenameVector, isRNA);
}
