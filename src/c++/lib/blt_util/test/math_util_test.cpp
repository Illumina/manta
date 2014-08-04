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

    BOOST_CHECK_CLOSE(result, expect, eps);
}


BOOST_AUTO_TEST_CASE( test_log_sum )
{

    single_log_sum_test(0.5,0.2);
    single_log_sum_test(0.00001,0.00000001);
    single_log_sum_test(1,1);
}

BOOST_AUTO_TEST_SUITE_END()

