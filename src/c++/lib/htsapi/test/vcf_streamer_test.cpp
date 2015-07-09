// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2015 Illumina, Inc.
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

#include "test_config.h"

#include "htsapi/vcf_streamer.hh"

#include "boost/test/unit_test.hpp"



BOOST_AUTO_TEST_SUITE( test_vcf_streamer )

static
const char*
getTestpath()
{
    static const std::string testPath(std::string(TEST_DATA_PATH) + "/vcf_streamer_test.vcf.gz");
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

