// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "boost/test/unit_test.hpp"

#include "htsapi/SimpleAlignment_bam_util.hh"
#include "manta/SVLocusScanner.cpp"
#include "truth/TruthTracker.hh"


BOOST_AUTO_TEST_SUITE( test_SVLocusScanner )


BOOST_AUTO_TEST_CASE( test_getSVCandidatesFromReadIndels )
{
    const ReadScannerOptions opt;
    const ReadScannerDerivOptions dopt(opt);

    ALIGNPATH::path_t inputPath;
    cigar_to_apath("100M2000D100M",inputPath);

    bam_record bamRead;
    bam1_t* bamDataPtr(bamRead.get_data());
    edit_bam_cigar(inputPath,*bamDataPtr);

    SimpleAlignment align(getAlignment(bamRead));

    std::vector<SVObservation> candidates;

    const std::string emptyPathStr;
    bam_header_info hdr_info;
    TruthTracker truthTracker(emptyPathStr, hdr_info);
    TrackedCandidates trackedCandidates(candidates,truthTracker);

    getSVCandidatesFromReadIndels(opt, dopt, align, trackedCandidates);

    BOOST_REQUIRE_EQUAL(candidates.size(),1u);
    BOOST_REQUIRE(candidates[0].bp1.interval.range.is_pos_intersect(100));
    BOOST_REQUIRE(candidates[0].bp2.interval.range.is_pos_intersect(2100));
}


BOOST_AUTO_TEST_SUITE_END()

