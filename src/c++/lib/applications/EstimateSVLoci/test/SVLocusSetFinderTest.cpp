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

#include "applications/EstimateSVLoci/SVLocusSetFinder.hpp"
#include "test/testAlignmentDataUtil.hpp"
#include "test/testFileMakers.hpp"

#include "boost/make_unique.hpp"

/// Utility function necessary for testing the SVLocusSetFinder.
/// Necessary to put here instead of EstimateSVLoci_util.hh because of undefined reference issue.
/// TODO refactor code to remove the necessity for this to exist.
static std::unique_ptr<SVLocusSetFinder> buildSVLocusSetFinder()
{
  const bam_header_info    bamHeaderInfo(buildTestBamHeader());
  TestStatsFileMaker       statsFile;
  TestAlignHeaderFileMaker alignFile(bamHeaderInfo);

  ESLOptions opts;
  opts.alignFileOpt.alignmentFilenames = {alignFile.getFilename()};
  opts.alignFileOpt.isAlignmentTumor   = {false};
  opts.statsFilename                   = statsFile.getFilename();

  return boost::make_unique<SVLocusSetFinder>(
      opts,
      GenomeInterval(0, 0, 499),
      std::make_shared<reference_contig_segment>(),
      std::make_shared<SVLocusSet>(opts.graphOpt, bamHeaderInfo, opts.alignFileOpt.alignmentFilenames));
}

BOOST_AUTO_TEST_SUITE(SVLocusSetFinderUpdate_test_suite)

#if 0
BOOST_AUTO_TEST_CASE( test_DepthFiltering )
{
    BOOST_TEST_MESSAGE("SDS MANTA-753");
    //TODO Setup a chrom depth for a chromosome and then setup read depth
}
#endif

BOOST_AUTO_TEST_CASE(test_MapQuality_Filtering)
{
  BOOST_TEST_MESSAGE("SDS MANTA-754");

  // TODO Setup read with quality above, below, and same as _opt.minMapq

  std::unique_ptr<SVLocusSetFinder> svLSF(buildSVLocusSetFinder());

  // stand-in state reporter used for the unit test, does nothing
  stream_state_reporter dummyStateReporter;

  // Valid anomalous read fails because its mapQ is 14 and minMapQ is 15.
  bam_record bamRead2;
  buildTestBamRecord(bamRead2, 0, 200, 0, 100, 99, 14);
  svLSF->update(dummyStateReporter, bamRead2, bamRead2.target_id());
  BOOST_REQUIRE_EQUAL(svLSF->getLocusSet().size(), 0);

  // Valid anomalous read passes because its mapQ is 15 and minMapQ is 15.
  bam_record bamRead1;
  buildTestBamRecord(bamRead1, 0, 200, 0, 100, 99, 15);
  svLSF->update(dummyStateReporter, bamRead1, bamRead1.target_id());
  BOOST_REQUIRE_EQUAL(svLSF->getLocusSet().size(), 1);

  svLSF->flush();

  BOOST_TEST_MESSAGE("SDS MANTA-755");

  // Test that the mapQ 14 is filtered and 15 is not filtered.
  //
  const auto& inputCounts(svLSF->getLocusSet().getAllSampleReadCounts().getSampleCounts(0).input);
  BOOST_REQUIRE_EQUAL(inputCounts.minMapq, 1);
  BOOST_REQUIRE_EQUAL(inputCounts.evidenceCount.total, 1);

  //// Indirect test matching above after write/read stats file:
  //const SVLocusSetStatsFileMaker stats(svLSF->getLocusSet());
  //BOOST_REQUIRE_EQUAL(getValueFromTSVKeyValFile(stats.getFilename(), "MinMapqFiltered"), "1");
  //BOOST_REQUIRE_EQUAL(getValueFromTSVKeyValFile(stats.getFilename(), "NotFiltered"), "1");
}

BOOST_AUTO_TEST_CASE(test_SplitReadSemiAligned)
{
  // create split read - not semi-aligned evidence
  bam_record splitRead;
  {
    static const pos_t alignPos(250);
    // Semi-aligned read sequence - is evidence
    static const char semiQuerySeq[] = "AACCCACAAACATCACACACAAGAGTCCAGAGACCGACTTTTTTCTAAAA";

    buildTestBamRecord(splitRead, 0, alignPos, 0, alignPos + 50, 8, 15, "50M", semiQuerySeq);
    addSupplementaryAlignmentEvidence(splitRead);
  }

  // insert split read in locus set finder:
  std::unique_ptr<SVLocusSetFinder> svLSF(buildSVLocusSetFinder());

  // stand-in state reporter used for the unit test, does nothing
  stream_state_reporter dummyStateReporter;
  svLSF->update(dummyStateReporter, splitRead, 0);
  svLSF->flush();

  // test locus set finder stats:
  // split read called first, not called for semi-aligned read sequence.
  const auto& inputCounts(svLSF->getLocusSet().getAllSampleReadCounts().getSampleCounts(0).input);
  BOOST_REQUIRE_EQUAL(inputCounts.evidenceCount.assm, 0);
  BOOST_REQUIRE_EQUAL(inputCounts.evidenceCount.split, 1);

  // repeats the same test as above after going through read/write cycle of stats reporting
  //const SVLocusSetStatsFileMaker stats(svLSF->getLocusSet());
  //BOOST_REQUIRE_EQUAL(getValueFromTSVKeyValFile(stats.getFilename(), "NotFilteredAndSemiAligned"), "0");
  //BOOST_REQUIRE_EQUAL(getValueFromTSVKeyValFile(stats.getFilename(), "NotFilteredAndSplitRead"), "1");
}

BOOST_AUTO_TEST_SUITE_END()
