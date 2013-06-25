// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

#include "boost/test/unit_test.hpp"

#include "bam_header_util.hh"


BOOST_AUTO_TEST_SUITE( test_bam_header_util )


BOOST_AUTO_TEST_CASE( test_parse_bam_region )
{

    bam_header_info header;

    header.chrom_data.push_back(bam_header_info::chrom_info("chr1",1000));
    header.chrom_data.push_back(bam_header_info::chrom_info("chr2",2000));

    int32_t tid,begin_pos,end_pos;
    parse_bam_region(header,std::string("chr2"),tid,begin_pos,end_pos);

    BOOST_REQUIRE_EQUAL(tid,1);
    BOOST_REQUIRE_EQUAL(begin_pos,0);
    BOOST_REQUIRE_EQUAL(end_pos,2000);

    parse_bam_region(header,std::string("chr1:100"),tid,begin_pos,end_pos);

    BOOST_REQUIRE_EQUAL(tid,0);
    BOOST_REQUIRE_EQUAL(begin_pos,99);
    BOOST_REQUIRE_EQUAL(end_pos,1000);

    parse_bam_region(header,std::string("chr1:100-200"),tid,begin_pos,end_pos);

    BOOST_REQUIRE_EQUAL(tid,0);
    BOOST_REQUIRE_EQUAL(begin_pos,99);
    BOOST_REQUIRE_EQUAL(end_pos,200);
}


BOOST_AUTO_TEST_SUITE_END()

