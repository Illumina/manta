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

/// \file
/// \author Trevor Ramsay
///

#include "boost/test/unit_test.hpp"

#include "BamStreamerUtils.hpp"

#include "test/testAlignmentDataUtil.hpp"
#include "test/testFileMakers.hpp"

#include "boost/filesystem.hpp"

BOOST_AUTO_TEST_SUITE(BamStreamerUtils_test_suite)

struct OpenBamStreamsFixture {
  OpenBamStreamsFixture() : readValues({{0, 10}, {0, 130}, {1, 20}, {1, 140}})
  {
    const bam_header_info   bamHeader(buildTestBamHeader());
    std::vector<bam_record> readsToAdd;
    for (const auto& readVal : readValues) {
      readsToAdd.emplace_back();
      buildTestBamRecord(readsToAdd.back(), readVal.first, readVal.second);
    }
    buildTestBamFile(bamHeader, readsToAdd, _bamFilename());

    const std::string        referenceFilename = getTestReferenceFilename();
    std::vector<std::string> bamFilenames      = {_bamFilename()};
    openBamStreams(referenceFilename, bamFilenames, bamStreams);
  }

  const std::vector<std::pair<int, int>>     readValues;
  std::vector<std::shared_ptr<bam_streamer>> bamStreams;

private:
  const std::string&     _bamFilename() const { return _bamFilenameMaker.getFilename(); }
  const BamFilenameMaker _bamFilenameMaker;
};

BOOST_FIXTURE_TEST_SUITE(openBamStreams_test_suite, OpenBamStreamsFixture)

BOOST_AUTO_TEST_CASE(test_openBamStreams)
{
  auto& bamStream(*(bamStreams[0]));
  for (const auto& readVal : readValues) {
    BOOST_REQUIRE(bamStream.next());
    const bam_record& read(*(bamStream.get_record_ptr()));
    BOOST_REQUIRE_EQUAL(read.pos() - 1, readVal.second);
  }
  BOOST_REQUIRE(!bamStream.next());
}

BOOST_AUTO_TEST_CASE(test_openBamStreamsRegion)
{
  resetBamStreamsRegion("chrFoo", bamStreams);

  auto& bamStream(*(bamStreams[0]));
  for (const auto& readVal : readValues) {
    if (readVal.first != 0) continue;
    BOOST_REQUIRE(bamStream.next());
    const bam_record& read(*(bamStream.get_record_ptr()));
    BOOST_REQUIRE_EQUAL(read.pos() - 1, readVal.second);
  }
  BOOST_REQUIRE(!bamStream.next());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
