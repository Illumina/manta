//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2017 Illumina, Inc.
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

///
/// \file
/// \author Trevor Ramsay
///

#include "boost/test/unit_test.hpp"

#include "manta/test/testUtil.hh"
#include "applications/EstimateSVLoci/SVLocusSetFinder.hh"


/// Utility function necessary for testing the SVLocusSetFinder.
/// Necessary to put here instead of EstimateSVLoci_util.hh because of undefined reference issue.
/// TODO refactor code to remove the necessity for this to exist.
static std::unique_ptr<SVLocusSetFinder> buildSVLocusSetFinder(
    const bam_header_info& bamHeaderInput,
    const std::string& alignFilename,
    const std::string& statsFilename = std::string("tempStatsFile.txt"),
    const GenomeInterval& interval = GenomeInterval(1, 1, 1000))
{
    ESLOptions opts;
    reference_contig_segment seg;

    createStatsFile(statsFilename);
    createAlignFile(bamHeaderInput, alignFilename);

    std::vector<bool> tempVector(1, false);

    std::vector<std::string> alignFilenameVector(1, alignFilename);
    AlignmentFileOptions afo;

    afo.alignmentFilename = alignFilenameVector;
    afo.isAlignmentTumor = tempVector;

    opts.alignFileOpt = afo;
    opts.statsFilename = statsFilename;

    std::unique_ptr<SVLocusSetFinder> newSVLSF(new SVLocusSetFinder(opts, interval, bamHeaderInput, seg));

    return newSVLSF;
}

BOOST_AUTO_TEST_SUITE( test_SVLocusSetFinderUpdate )

BOOST_AUTO_TEST_CASE( test_DepthFiltering )
{
    //BOOST_TEST_MESSAGE("SDS MANTA-753");
    //TODO Setup a chrom depth for a chromosome and then setup read depth
}

BOOST_AUTO_TEST_CASE( test_MapQuality_Filtering )
{
    BOOST_TEST_MESSAGE("SDS MANTA-754");
    //TODO Setup read with quality above, below, and same as _opt.minMapq
    const std::string& alignFilename = std::string("tempAlignFile.txt");
    bam_header_info bamHeader = bam_header_info();

    bamHeader.chrom_to_index.insert(std::pair<std::string,int32_t>("chrM",1000000));
    bamHeader.chrom_data.emplace_back("chrM",1000000);

    std::unique_ptr<SVLocusSetFinder> svLSF(buildSVLocusSetFinder(bamHeader, alignFilename));

    // stand-in state reporter used for the unit test, does nothing
    stream_state_reporter dummyStateReporter;

    // Valid anomolous read fails because its mapQ is 14 and minMapQ is 15.
    bam_record bamRead2;
    buildBamRecord(bamRead2, 0, 200, 0, 100, 99, 14);
    svLSF->update(dummyStateReporter, bamRead2, bamRead2.target_id());
    BOOST_REQUIRE_EQUAL(svLSF->getLocusSet().size(), 0);

    // Valid anomolous read passes because its mapQ is 15 and minMapQ is 15.
    bam_record bamRead1;
    buildBamRecord(bamRead1, 0, 200, 0, 100, 99, 15);
    svLSF->update(dummyStateReporter, bamRead1, bamRead1.target_id());
    BOOST_REQUIRE_EQUAL(svLSF->getLocusSet().size(), 1);

    BOOST_TEST_MESSAGE("SDS MANTA-755");
    std::string statsFile = "dumpStatsFile.txt";
    std::filebuf fb;
    fb.open(statsFile, std::ios::out);
    std::ostream os(&fb);

    svLSF->getLocusSet().dumpStats(os);
    fb.close();

    // Test that the mapQ 14 is filtered and 15 is not filtered.
    BOOST_REQUIRE_EQUAL(getResultFromStatsFile(statsFile, "MinMapqFiltered"), "1");
    BOOST_REQUIRE_EQUAL(getResultFromStatsFile(statsFile, "NotFiltered"), "1");

    remove(alignFilename.c_str());
    remove(statsFile.c_str());
}

BOOST_AUTO_TEST_CASE( test_SplitReadSemiAligned )
{
    static const pos_t alignPos(500);

    // minCandidateVariantSize = 8
    bam_header_info bamHeader = buildBamHeader();
    std::unique_ptr<SVLocusScanner> scanner(buildSVLocusScanner(bamHeader));

    //const unsigned defaultReadGroupIndex = 0;
    reference_contig_segment ref = reference_contig_segment();
    static const char refSeq[]   =      "AACCTTTTTTCATCACACACAAGAGTCCAGAGACCGACTTCCCCCCAAAA";
    ref.seq() = refSeq;
    ref.set_offset(alignPos);
    //SVLocusEvidenceCount *incountsPtr = new SVLocusEvidenceCount();

    // Semi-aligned read sequence - is evidence
    static const char semiQuerySeq[]   = "AACCCACAAACATCACACACAAGAGTCCAGAGACCGACTTTTTTCTAAAA";

    // split read - not semi-aligned evidence
    std::string dumpStatsFile = "dumpStatsFile.txt";
    bam_record splitRead;
    buildBamRecord(splitRead, 0, alignPos, 0, alignPos+50, 8, 15, "50M", semiQuerySeq);
    addSupplementaryEvidence(splitRead);

    // stand-in state reporter used for the unit test, does nothing
    stream_state_reporter dummyStateReporter;

    const std::string& alignFilename = std::string("tempAlignFile.txt");
    std::unique_ptr<SVLocusSetFinder> finder(buildSVLocusSetFinder(bamHeader, alignFilename));
    finder->update(dummyStateReporter, splitRead, 0);

    std::filebuf fb;
    fb.open(dumpStatsFile, std::ios::out);
    std::ostream os(&fb);

    finder->getLocusSet().dumpStats(os);
    fb.close();
    // split read called first, not called for semi-aligned read sequence.
    BOOST_REQUIRE_EQUAL(getResultFromStatsFile(dumpStatsFile, "NotFilteredAndSemiAligned"), "0");
    BOOST_REQUIRE_EQUAL(getResultFromStatsFile(dumpStatsFile, "NotFilteredAndSplitRead"), "1");

    remove(dumpStatsFile.c_str());
}

BOOST_AUTO_TEST_SUITE_END()
