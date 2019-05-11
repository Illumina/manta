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
/// \author Trevor Ramsay
///

#include "boost/test/unit_test.hpp"

#include "EstimateSVLociRunner.hpp"

#include "test/testAlignmentDataUtil.hpp"
#include "test/testFileMakers.hpp"
#include "test/testUtil.hpp"

#include "boost/make_unique.hpp"

/// \brief Construct an EstimateSVLociRunner with dummy options and small dummy data files
struct ConstructTestEstimateSVLociRunner {
  ConstructTestEstimateSVLociRunner()
  {
    // Programatically construct a bam file for testing
    std::vector<std::string> bamFiles = {_bamFilename.getFilename()};
    {
      std::vector<bam_record> readsToAdd(2);

      // Valid anomalous read passes because its mapQ is 15 and minMapQ is 15.
      buildTestBamRecord(readsToAdd[0], 0, 200, 0, 100, 100, 15);

      // Valid anomalous read fails because its mapQ is 14 and minMapQ is 15.
      buildTestBamRecord(readsToAdd[1], 0, 300, 0, 150, 100, 14);

      buildTestBamFile(buildTestBamHeader(), readsToAdd, bamFiles[0]);
    }

    // Initialize dummy ESL options:
    ESLOptions eslOpt;
    eslOpt.referenceFilename = getTestReferenceFilename();
    eslOpt.statsFilename     = _statsFileMaker.getFilename();

    eslOpt.alignFileOpt.alignmentFilenames = bamFiles;
    eslOpt.alignFileOpt.isAlignmentTumor   = {false};

    _eslRunnPtr = boost::make_unique<EstimateSVLociRunner>(eslOpt);
  }

  EstimateSVLociRunner& get() { return *_eslRunnPtr; }

private:
  BamFilenameMaker                      _bamFilename;
  TestStatsFileMaker                    _statsFileMaker;
  std::unique_ptr<EstimateSVLociRunner> _eslRunnPtr;
};

BOOST_AUTO_TEST_SUITE(EstimateSVLociRunner_test_suite)

// Test EstimateSVLociRunner's construction of the full SVLocusSet object in its constructor.
BOOST_AUTO_TEST_CASE(test_SVLocusSetInitialization)
{
  ConstructTestEstimateSVLociRunner testESLRunner;
  auto&                             eslRunner(testESLRunner.get());

  BOOST_REQUIRE_EQUAL(eslRunner.getLocusSet().getAllSampleReadCounts().size(), 1);
}

// Test that the estimateSVLociForSingleRegion works correctly with a basic input of a read filtered for mapQ
// and a read that should not be filtered at all. The results will be output to a stats file for confirmation
// to satisfy MANTA-755
BOOST_AUTO_TEST_CASE(test_SVLocusSampleCounts)
{
  BOOST_TEST_MESSAGE("SDS MANTA-755");

  ConstructTestEstimateSVLociRunner testESLRunner;
  auto&                             eslRunner(testESLRunner.get());

  eslRunner.estimateSVLociForSingleRegion("chrFoo");

  // Test that the mapQ 14 is filtered and 15 is not filtered.
  // 1. direct test:
  const auto& inputCounts(eslRunner.getLocusSet().getAllSampleReadCounts().getSampleCounts(0).input);
  BOOST_REQUIRE_EQUAL(inputCounts.minMapq, 1);
  BOOST_REQUIRE_EQUAL(inputCounts.evidenceCount.total, 1);

  // 2. indirect test through stats file
  SVLocusSetStatsFileMaker graphStats(eslRunner.getLocusSet());
  BOOST_REQUIRE_EQUAL(getValueFromTSVKeyValFile(graphStats.getFilename(), "MinMapqFiltered"), "1");
  BOOST_REQUIRE_EQUAL(getValueFromTSVKeyValFile(graphStats.getFilename(), "NotFiltered"), "1");
}

BOOST_AUTO_TEST_SUITE_END()
