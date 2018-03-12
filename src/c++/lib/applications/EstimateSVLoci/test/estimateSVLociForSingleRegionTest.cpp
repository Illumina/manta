//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2018 Illumina, Inc.
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

#include "estimateSVLociForSingleRegion.hh"

#include "test/testAlignmentDataUtil.hh"
#include "test/testFileMakers.hh"

#include <fstream>

BOOST_AUTO_TEST_SUITE( estimateSVLociForSingleRegion_test_suite )

/// Test that the estimateSVLociForSingleRegion works correctly with a basic input of a read filtered for mapQ and a
/// read that should not be filtered at all. The results will be output to a stats file for confirmation to satisfy
/// MANTA-755
BOOST_AUTO_TEST_CASE( test_SVLocusSampleCounts )
{
    BOOST_TEST_MESSAGE("SDS MANTA-755");

    // Programatically construct a bam file for testing
    TestFilenameMaker bamFilename;
    std::vector<std::string> bamFiles = { bamFilename.getFilename() };
    {
        std::vector<bam_record> readsToAdd(2);

        // Valid anomalous read passes because its mapQ is 15 and minMapQ is 15.
        buildTestBamRecord(readsToAdd[0], 0, 200, 0, 100, 100, 15);

        // Valid anomalous read fails because its mapQ is 14 and minMapQ is 15.
        buildTestBamRecord(readsToAdd[1], 0, 300, 0, 150, 100, 14);

        buildTestBamFile(buildTestBamHeader(), readsToAdd, bamFiles[0]);
    }

    // Initialize estimateSVLoci with dummy options

    TestStatsFileMaker statsFileMaker;

    // Initialize fake ESL options:
    ESLOptions eslOpt;
    eslOpt.referenceFilename = getTestReferenceFilename();
    eslOpt.statsFilename = statsFileMaker.getFilename();

    eslOpt.alignFileOpt.alignmentFilenames = bamFiles;
    eslOpt.alignFileOpt.isAlignmentTumor = { false };

    SVLocusSet svLoci;
    estimateSVLociForSingleRegion(eslOpt, "chrFoo", svLoci);

    // Test that the mapQ 14 is filtered and 15 is not filtered.
    //
    const auto& sampleCounts(svLoci.getCounts().getSampleCounts(0));
    BOOST_REQUIRE_EQUAL(sampleCounts.input.minMapq, 1);
    BOOST_REQUIRE_EQUAL(sampleCounts.input.evidenceCount.total, 1);

    // Test that the mapQ 14 is filtered and 15 is not filtered.
    //SVLocusSetStatsFileMaker graphStats(svLoci);
    //BOOST_REQUIRE_EQUAL(getValueFromTSVKeyValFile(graphStats.getFilename(), "MinMapqFiltered"), "1");
    //BOOST_REQUIRE_EQUAL(getValueFromTSVKeyValFile(graphStats.getFilename(), "NotFiltered"), "1");
}

BOOST_AUTO_TEST_SUITE_END()
