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

#include "boost/test/unit_test.hpp"

#include "depth_buffer.hh"


BOOST_AUTO_TEST_SUITE( test_depth_buffer )


/// return buffer loaded with simple test pattern
///
/// at depth 10Y, depth is Y
///
static
depth_buffer
get_db_test_pattern()
{
    depth_buffer db;

    // load a depth pattern in:
    for (unsigned i(100); i<110; ++i)
    {
        for (unsigned j(i+1); j<110; ++j)
        {
            db.inc(j);
        }
    }
    return db;
}


/// return buffer loaded with simple test pattern
///
/// at pos Y, depth is Y-100 before compression
///
static
depth_buffer_compressible
get_db_compressible_test_pattern(
    const unsigned compressionLevel)
{
    depth_buffer_compressible db(compressionLevel);

    // load a depth pattern in:
    for (unsigned i(101); i<200; ++i)
    {
        db.inc(i,200-i);
    }
    return db;
}


BOOST_AUTO_TEST_CASE( test_depth_buffer_val )
{
    depth_buffer db(get_db_test_pattern());
    BOOST_CHECK_EQUAL(static_cast<int>(db.val(109)),9);
}


BOOST_AUTO_TEST_CASE( test_depth_buffer_clear )
{
    depth_buffer db(get_db_test_pattern());
    db.clear_pos(109);
    BOOST_CHECK_EQUAL(static_cast<int>(db.val(109)),0);
}


BOOST_AUTO_TEST_CASE( test_depth_buffer_range )
{
    depth_buffer db(get_db_test_pattern());
    BOOST_CHECK(! db.is_range_ge_than(0,107,8));
    BOOST_CHECK(  db.is_range_ge_than(0,108,8));
}

BOOST_AUTO_TEST_CASE( test_depth_buffer_compressible_val )
{
    depth_buffer_compressible db(get_db_compressible_test_pattern(8));
    BOOST_CHECK_EQUAL(static_cast<int>(db.val(107)),8);
    BOOST_CHECK_EQUAL(static_cast<int>(db.val(110)),8);
    BOOST_CHECK_EQUAL(static_cast<int>(db.val(150)),48);
    BOOST_CHECK_EQUAL(static_cast<int>(db.val(199)),96);
}


BOOST_AUTO_TEST_CASE( test_depth_buffer_compressible_clear )
{
    depth_buffer_compressible db(get_db_compressible_test_pattern(8));
    for (unsigned i(100); i<120; ++i)
    {
        db.clear_pos(i);
    }
    BOOST_CHECK_EQUAL(static_cast<int>(db.val(110)),0);
    BOOST_CHECK_EQUAL(static_cast<int>(db.val(150)),48);
    BOOST_CHECK_EQUAL(static_cast<int>(db.val(199)),96);
    BOOST_CHECK_EQUAL(static_cast<int>(db.val(109)),0);
}


BOOST_AUTO_TEST_SUITE_END()

