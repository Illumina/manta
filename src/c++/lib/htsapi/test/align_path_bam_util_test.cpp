// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

#include "htsapi/align_path_bam_util.hh"
#include "htsapi/bam_record.hh"

#include "boost/test/unit_test.hpp"



BOOST_AUTO_TEST_SUITE( test_align_path_bam_util )


BOOST_AUTO_TEST_CASE( test_edit_bam_cigar )
{
    const std::string testCigar("10M10D1I10M");
    ALIGNPATH::path_t inputPath;
    cigar_to_apath(testCigar.c_str(),inputPath);

    bam_record bamRead;
    bam1_t* bamDataPtr(bamRead.get_data());
    edit_bam_cigar(inputPath,*bamDataPtr);

    ALIGNPATH::path_t outputPath;
    bam_cigar_to_apath(bamRead.raw_cigar(),bamRead.n_cigar(),outputPath);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(outputPath),testCigar);
}


BOOST_AUTO_TEST_SUITE_END()

