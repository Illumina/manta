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

#include <iostream>
#include "boost/test/unit_test.hpp"

#include "SVEvidence.hpp"

BOOST_AUTO_TEST_SUITE(SVEvidence_test_suite)

// Comparing all fields of SVFragmentEvidenceAlleleBreakendPerRead of observed data with expected
// SVFragmentEvidenceAlleleBreakendPerRead data. Following fields are compared:
// 1. isSplitSupport - Does this read support this allele in the breakpoint?
// 2. isSplitEvaluated -  have we checked this read for split support of this breakpoint?
// 3. isTier2SplitSupport - Does this read support this allele in the breakpoint by permissive criteria?
// 4. splitEvidence - Evidence score
// 5. splitLnLhood - ln likelihood of best split alignment
static void compareAlleleBreakend(
    SVFragmentEvidenceAlleleBreakendPerRead observed, SVFragmentEvidenceAlleleBreakendPerRead expected)
{
  BOOST_REQUIRE_EQUAL(observed.isSplitSupport, expected.isSplitSupport);
  BOOST_REQUIRE_EQUAL(observed.isSplitEvaluated, expected.isSplitEvaluated);
  BOOST_REQUIRE_EQUAL(observed.isTier2SplitSupport, expected.isTier2SplitSupport);
  BOOST_REQUIRE_EQUAL(observed.splitEvidence, expected.splitEvidence);
  BOOST_REQUIRE_EQUAL(observed.splitLnLhood, expected.splitLnLhood);
}

// Comparing all fields of SVFragmentEvidenceRead of observed data with expected
// SVFragmentEvidenceRead data. Following fields are compared:
// 1. isScanned - Whether this read is scanned or not.
// 2. isShadow - read was originally unmapped but had a mapped mate read, mapq is MAPQ of the mate in this
// case
// 3. mapq - mapping quality
// 4. size - Size of the data
static void compareSVEvidenceRead(SVFragmentEvidenceRead observed, SVFragmentEvidenceRead expected)
{
  BOOST_REQUIRE_EQUAL(observed.isScanned, expected.isScanned);
  BOOST_REQUIRE_EQUAL(observed.isShadow, expected.isShadow);
  BOOST_REQUIRE_EQUAL(observed.mapq, expected.mapq);
  BOOST_REQUIRE_EQUAL(observed.size, expected.size);
}

// Test the SVFragmentEvidenceAlleleBreaked structure
// Following points need to be tested
// 1. Whether value of a field in the object matches with the value returned by the api.
// 2. If there is any clear() api, it should reset the value of a variable to default.
//
BOOST_AUTO_TEST_CASE(test_SVFragmentEvidenceAlleleBreakend)
{
  SVFragmentEvidenceAlleleBreakend fragmentEvidenceAltBp1;

  // Setting the values for read1 of svFragmentEvidenceAlleleBreakend
  fragmentEvidenceAltBp1.read1.isSplitEvaluated    = true;
  fragmentEvidenceAltBp1.read1.isSplitSupport      = true;
  fragmentEvidenceAltBp1.read1.isTier2SplitSupport = true;
  fragmentEvidenceAltBp1.read1.splitEvidence       = 0.5;
  fragmentEvidenceAltBp1.read1.splitLnLhood        = 0.2;

  // setting the values for read2 of svFragmentEvidenceAlleleBreakend
  fragmentEvidenceAltBp1.read2.isSplitEvaluated    = true;
  fragmentEvidenceAltBp1.read2.isSplitSupport      = false;
  fragmentEvidenceAltBp1.read2.isTier2SplitSupport = false;
  fragmentEvidenceAltBp1.read2.splitEvidence       = 0.1;
  fragmentEvidenceAltBp1.read2.splitLnLhood        = 0.05;

  // Check the attributes of read1 and read2 (returned from api) whether they are
  // matching with the required values
  compareAlleleBreakend(fragmentEvidenceAltBp1.getRead(true), fragmentEvidenceAltBp1.read1);
  compareAlleleBreakend(fragmentEvidenceAltBp1.getRead(false), fragmentEvidenceAltBp1.read2);

  // Setting the value of fragmentLengthProb which is probability of the fragment size given
  // this allele. Also setting the isFragmentSupport which is whether this read pair fragment supports
  // this allele on this breakend or not. Checking whether values are set correctly or not and after
  // that calling clearPairSupport to check whether the state goes back to default or not.
  fragmentEvidenceAltBp1.fragLengthProb    = 0.5;
  fragmentEvidenceAltBp1.isFragmentSupport = true;
  BOOST_REQUIRE(fragmentEvidenceAltBp1.isFragmentSupport);
  BOOST_REQUIRE_EQUAL(fragmentEvidenceAltBp1.fragLengthProb, 0.5);

  // clearing means it should reset to default values.
  fragmentEvidenceAltBp1.clearPairSupport();
  BOOST_REQUIRE(!fragmentEvidenceAltBp1.isFragmentSupport);
  BOOST_REQUIRE_EQUAL(fragmentEvidenceAltBp1.fragLengthProb, 0);

  // Setting the values for read1 of SVFragmentEvidenceRead
  SVFragmentEvidenceRead read1;
  read1.isScanned = true;
  read1.isShadow  = true;
  read1.size      = 150;
  read1.mapq      = 70;
  // Setting the values for read2 of SVFragmentEvidenceRead
  SVFragmentEvidenceRead read2;
  read2.isShadow  = false;
  read2.isScanned = false;
  read2.size      = 100;
  read2.mapq      = 40;

  // Check the attributes of read1 and read2 (returned from api) whether they are matching with
  // the required values
  SVFragmentEvidence svFragmentEvidence;
  svFragmentEvidence.read1 = read1;
  svFragmentEvidence.read2 = read2;
  compareSVEvidenceRead(svFragmentEvidence.getRead(true), read1);
  compareSVEvidenceRead(svFragmentEvidence.getRead(false), read2);

  fragmentEvidenceAltBp1.fragLengthProb    = 0.5;
  fragmentEvidenceAltBp1.isFragmentSupport = true;
  svFragmentEvidence.alt.bp1               = fragmentEvidenceAltBp1;

  // Whether this fragment read provides any pair evidence for any breakpoint of the ALT allele or not.
  // That means isFragmentSupport of either bp1 or bp2 should be true. Here
  // isFragmentSupport of bp1 is true. So this should return true.
  BOOST_REQUIRE(svFragmentEvidence.isAltSpanningPairSupport());
  svFragmentEvidence.alt.bp1.isFragmentSupport = false;
  svFragmentEvidence.alt.bp2.isFragmentSupport = true;
  // Here isFragmentSupport of bp2 is true but bp1 is false. So this should return true.
  BOOST_REQUIRE(svFragmentEvidence.isAltSpanningPairSupport());
  svFragmentEvidence.alt.bp1.isFragmentSupport = true;
  svFragmentEvidence.alt.bp2.isFragmentSupport = true;
  // Here isFragmentSupport of both bp1 and bp2 are true. So this should return true.
  BOOST_REQUIRE(svFragmentEvidence.isAltSpanningPairSupport());
  svFragmentEvidence.alt.bp1.isFragmentSupport = false;
  svFragmentEvidence.alt.bp2.isFragmentSupport = false;
  // Here isFragmentSupport of both bp1 and bp2 are false. So this should return false.
  BOOST_REQUIRE(!svFragmentEvidence.isAltSpanningPairSupport());

  // Whether this fragment provides any pair evidence for any allele/ref combination or not.
  // That means this method will return true if either isAltSpanningPairSupport is true
  // or ref svEvidence is true.
  // ref svEvidence should be true if isFragmentSupport of either bp1 or bp2 is true.
  // Here isFragmentSupport of both bp1 and bp2 are true. So this should return true.
  svFragmentEvidence.ref.bp1.isFragmentSupport = true;
  svFragmentEvidence.ref.bp2.isFragmentSupport = true;
  BOOST_REQUIRE(svFragmentEvidence.isAnySpanningPairSupport());
  // Here isFragmentSupport of bp1 is true but bp2 is false. So this should return true.
  svFragmentEvidence.ref.bp1.isFragmentSupport = true;
  svFragmentEvidence.ref.bp2.isFragmentSupport = false;
  BOOST_REQUIRE(svFragmentEvidence.isAnySpanningPairSupport());
  // Here isFragmentSupport of bp2 is true but bp1 is false. So this should return true.
  svFragmentEvidence.ref.bp1.isFragmentSupport = false;
  svFragmentEvidence.ref.bp2.isFragmentSupport = true;
  BOOST_REQUIRE(svFragmentEvidence.isAnySpanningPairSupport());
  // Here isFragmentSupport of both bp1 and bp2 are false. So this should return false.
  svFragmentEvidence.ref.bp1.isFragmentSupport = false;
  svFragmentEvidence.ref.bp2.isFragmentSupport = false;
  BOOST_REQUIRE(!svFragmentEvidence.isAnySpanningPairSupport());

  // Whether this fragment read provides any split evidence for any bp of the ALT allele or not.
  // So isAltSplitReadSupport of read1 is true if and only if any of the breakpoint
  // supports split evidence.. Here breakpoint-1 of read-1 supports split evidence.
  BOOST_REQUIRE(svFragmentEvidence.isAltSplitReadSupport(true));
  BOOST_REQUIRE(svFragmentEvidence.isAnySplitReadSupport(true).first);
  BOOST_REQUIRE(!svFragmentEvidence.isAnySplitReadSupport(true).second);

  // Here breakpoint-1 and breakpoint-2 of read-2 do not support split evidence.
  BOOST_REQUIRE(!svFragmentEvidence.isAltSplitReadSupport(false));
  BOOST_REQUIRE(!svFragmentEvidence.isAnySplitReadSupport(false).first);
  BOOST_REQUIRE(!svFragmentEvidence.isAnySplitReadSupport(false).second);

  // Whether this fragment read supports any allele in the breakpoint by permissive criteria or not.
  // So isAltTier2SplitReadSupport of read1 is true if and only if any of the breakpoint
  // supports tier2 split. Here breakpoint-1 of read-1 supports tier2 split evidence.
  BOOST_REQUIRE(svFragmentEvidence.isAltTier2SplitReadSupport(true));
  BOOST_REQUIRE(svFragmentEvidence.isAnyTier2SplitReadSupport(true).first);
  BOOST_REQUIRE(!svFragmentEvidence.isAnyTier2SplitReadSupport(true).second);

  // Here breakpoint-1 and breakpoint-2 of read-2 do not support tier2 split.
  BOOST_REQUIRE(!svFragmentEvidence.isAltTier2SplitReadSupport(false));
  BOOST_REQUIRE(!svFragmentEvidence.isAnyTier2SplitReadSupport(false).first);
  BOOST_REQUIRE(!svFragmentEvidence.isAnyTier2SplitReadSupport(false).second);
}

BOOST_AUTO_TEST_SUITE_END()
