// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

#include "boost/test/unit_test.hpp"

#include "math_util.hh"

#include <cmath>


BOOST_AUTO_TEST_SUITE( math_util )

static
void
single_log_sum_test(const double x1,
                    const double x2)
{
    static const double eps(0.00001);

    const double expect(std::log(x1+x2));

    const double lnx1(std::log(x1));
    const double lnx2(std::log(x2));
    const double result(log_sum(lnx1,lnx2));

    BOOST_REQUIRE_CLOSE(result, expect, eps);
}


BOOST_AUTO_TEST_CASE( test_log_sum )
{

    single_log_sum_test(0.5,0.2);
    single_log_sum_test(0.00001,0.00000001);
    single_log_sum_test(1,1);
}


static
void
medianChecker(
    const int trueVal,
    const std::vector<int>& nums)
{
    std::vector<int> copy1(nums);
    BOOST_REQUIRE_EQUAL(trueVal,_ps_median(std::begin(copy1), std::end(copy1)));
    std::vector<int> copy2(nums);
    BOOST_REQUIRE_EQUAL(trueVal,_ne_median(std::begin(copy2), std::end(copy2)));
}


BOOST_AUTO_TEST_CASE( test_median )
{
    {
        std::vector<int> nums = {9,8,7,0,2,1,3};
        medianChecker(3,nums);
    }
    {
        std::vector<int> nums  = {9,8,7,0,2,1,3,10,11};
        medianChecker(7,nums);
    }
}


BOOST_AUTO_TEST_SUITE_END()

