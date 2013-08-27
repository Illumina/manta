// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
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
get_db_test_pattern() {
    depth_buffer db;

    // load a depth pattern in:
    for (unsigned i(100); i<110; ++i) {
        for (unsigned j(i+1); j<110; ++j) {
            db.inc(j);
        }
    }
    return db;
}


BOOST_AUTO_TEST_CASE( test_depth_buffer_val ) {

    depth_buffer db(get_db_test_pattern());
    BOOST_CHECK_EQUAL(static_cast<int>(db.val(109)),9);
}


BOOST_AUTO_TEST_CASE( test_depth_buffer_clear ) {

    depth_buffer db(get_db_test_pattern());
    db.clear_pos(109);
    BOOST_CHECK_EQUAL(static_cast<int>(db.val(109)),0);
}


BOOST_AUTO_TEST_CASE( test_depth_buffer_range ) {

    depth_buffer db(get_db_test_pattern());
    BOOST_CHECK(! db.is_range_ge_than(0,107,8));
    BOOST_CHECK(  db.is_range_ge_than(0,108,8));
}


BOOST_AUTO_TEST_SUITE_END()

