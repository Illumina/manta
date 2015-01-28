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

#include "boost/test/unit_test.hpp"

#include "io_util.hh"

#include <iomanip>
#include <sstream>


BOOST_AUTO_TEST_SUITE( test_io_util )


BOOST_AUTO_TEST_CASE( test_StreamScoper )
{
    static const double f(0.123456789);
    static const std::string strf4("0.1235");
    static const std::string strf2("0.12");

    std::ostringstream oss;

    oss <<  std::fixed << std::setprecision(4);

    oss << f;
    BOOST_REQUIRE_EQUAL(strf4,oss.str());

    oss.str("");

    {
        StreamScoper ss(oss);
        oss <<  std::fixed << std::setprecision(2);
        oss << f;
        BOOST_REQUIRE_EQUAL(strf2,oss.str());

        oss.str("");
    }

    oss << f;
    BOOST_REQUIRE_EQUAL(strf4,oss.str());
}


BOOST_AUTO_TEST_SUITE_END()
