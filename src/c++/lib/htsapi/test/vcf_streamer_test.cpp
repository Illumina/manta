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

#include "htsapi/vcf_streamer.hh"

#include "boost/test/unit_test.hpp"



BOOST_AUTO_TEST_SUITE( test_vcf_streamer )

static
const char*
getTestpath()
{
    static const std::string testPath("@CMAKE_CURRENT_SOURCE_DIR@/vcf_streamer_test.vcf.gz");
    return testPath.c_str();
}


BOOST_AUTO_TEST_CASE( test_vcf_streamer_region )
{
    vcf_streamer vcfs(getTestpath(),"chr1:750000-822000");

    const vcf_record* vptr(nullptr);

    BOOST_REQUIRE( vcfs.next() );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    BOOST_REQUIRE( vptr->is_valid() );
    BOOST_REQUIRE( vptr->is_indel() );
    BOOST_REQUIRE( ! vptr->is_snv() );

    BOOST_REQUIRE_EQUAL(vptr->pos, 757807);
    BOOST_REQUIRE_EQUAL(vptr->ref,"CCCTGGCCAGCAGATCCACCCTGTCTATACTACCTG");

    BOOST_REQUIRE( vcfs.next() );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    BOOST_REQUIRE( vptr->is_valid() );
    BOOST_REQUIRE( ! vptr->is_indel() );
    BOOST_REQUIRE( vptr->is_snv() );
    BOOST_REQUIRE_EQUAL(vptr->pos, 758807);
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),2u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"T");

    BOOST_REQUIRE( vcfs.next() );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    BOOST_REQUIRE( vptr->is_valid() );
    BOOST_REQUIRE( vptr->is_indel() );
    BOOST_REQUIRE( ! vptr->is_snv() );
    BOOST_REQUIRE_EQUAL(vptr->pos, 821604);
    BOOST_REQUIRE_EQUAL(vptr->alt.size(),1u);
    BOOST_REQUIRE_EQUAL(vptr->alt[0],"TGCCCTTTGGCAGAGCAGGTGTGCTGTGCTG");

    BOOST_REQUIRE( ! vcfs.next() );
    vptr = vcfs.get_record_ptr();
    assert(vptr == nullptr);
}

#if 0
BOOST_AUTO_TEST_CASE( test_vcf_streamer_noregion )
{
    vcf_streamer vcfs(getTestpath());

    const vcf_record* vptr(nullptr);

    BOOST_REQUIRE( vcfs.next() );
    vptr = vcfs.get_record_ptr();
    assert(vptr != nullptr);

    BOOST_REQUIRE( vptr->is_valid() );
    BOOST_REQUIRE( vptr->is_indel() );
    BOOST_REQUIRE( ! vptr->is_snv() );

    BOOST_REQUIRE_EQUAL(vptr->pos, 54712);

    BOOST_REQUIRE( vcfs.next() );
}
#endif

BOOST_AUTO_TEST_SUITE_END()

