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

#include "compat_util.hh"

#include <string>


BOOST_AUTO_TEST_SUITE( compat_util )

static
void
single_test_round(const double input,
                  const double expect)
{

    static const double eps(0.00001);
    BOOST_CHECK_CLOSE(compat_round(input), expect, eps);
}


BOOST_AUTO_TEST_CASE( test_round )
{
    single_test_round(3.5,4.0);
    single_test_round(3.2,3.0);
    single_test_round(3.7,4.0);
    single_test_round(-1.0,-1.0);
    single_test_round(-1.2,-1.0);
    single_test_round(-1.5,-2.0);
    single_test_round(-1.7,-2.0);
}



static
void
single_test_basename(const char* input,
                     const char* expect)
{

    const char* result(compat_basename(input));
    BOOST_CHECK_EQUAL(std::string(result), std::string(expect));
}


BOOST_AUTO_TEST_CASE( test_basename )
{
    single_test_basename("foo","foo");
    single_test_basename("","");

#ifdef _WIN32
    single_test_basename("\\foo","foo");
    single_test_basename("\\\\","");
    single_test_basename("\\","");
#else
    single_test_basename("/foo","foo");
    single_test_basename("//","");
    single_test_basename("/","");
#endif
}


BOOST_AUTO_TEST_SUITE_END()

