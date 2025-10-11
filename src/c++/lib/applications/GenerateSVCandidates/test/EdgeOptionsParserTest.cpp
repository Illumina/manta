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

#include "boost/test/unit_test.hpp"

#include "EdgeOptionsParser.hpp"

using namespace std;

BOOST_AUTO_TEST_SUITE(EdgeOptionsParser_test_suite)

// Test the parsing of edge options.
// Following points need to be tested
// 1. locus index format
// 2. bin count and bin index validation
BOOST_AUTO_TEST_CASE(test_ParseEdgeOptions)
{
  using namespace boost::program_options;

  EdgeOptions                           edgeoption;
  const char                            locusIndexKey[] = "locus-index";
  std::string                           errorMsg;
  boost::program_options::variables_map vm;
  string                                locusIndex = "1:2:3:4";
  vm.insert(std::make_pair(locusIndexKey, variable_value(locusIndex, true)));

  // locus index format is not correct. It should not have more than 3 colon separated segments.
  // locus index format is either of the following three types:
  // 1. locusIndex (single locus index)
  // 2. locusIndex:nodeIndex
  // 3. locusIndex:nodeIndex:nodeIndex
  parseOptions(vm, edgeoption, errorMsg);
  BOOST_REQUIRE(errorMsg != "" && errorMsg.find("locus-index") == 0);

  locusIndex = "1:2:3";  // locus-index format as locusIndex:nodeIndex:nodeIndex
  vm.clear();
  vm.insert(std::make_pair(locusIndexKey, variable_value(locusIndex, true)));
  edgeoption.binCount = 0;

  // Divide all edges in the graph into binCount bins of approx equal complexity.
  // So that at least 1 bin count is required.
  parseOptions(vm, edgeoption, errorMsg);
  BOOST_REQUIRE(errorMsg != "" && errorMsg.find("bin-count") == 0);

  edgeoption.binCount = 1;
  edgeoption.binIndex = 1;
  // Out of binCount bins, iterate through the edges in this bin only.
  // So the bin index should be less than total bin count.
  parseOptions(vm, edgeoption, errorMsg);
  BOOST_REQUIRE(errorMsg != "" && errorMsg.find("bin-index") == 0);

  // Successful case when all parameters are correct
  edgeoption.binCount = 1;
  edgeoption.binIndex = 0;
  parseOptions(vm, edgeoption, errorMsg);
  BOOST_REQUIRE_EQUAL(errorMsg, "");

  // Another successful case when locus index format is different.
  locusIndex = "1:2";  // locus-index format as locusIndex:nodeIndex
  vm.clear();
  vm.insert(std::make_pair(locusIndexKey, variable_value(locusIndex, true)));
  parseOptions(vm, edgeoption, errorMsg);
  BOOST_REQUIRE_EQUAL(errorMsg, "");

  // Another successful case when locus index format is different.
  locusIndex = "1";  // locus-index format as single locusIndex
  vm.clear();
  vm.insert(std::make_pair(locusIndexKey, variable_value(locusIndex, true)));
  parseOptions(vm, edgeoption, errorMsg);
  BOOST_REQUIRE_EQUAL(errorMsg, "");
}

BOOST_AUTO_TEST_SUITE_END()
