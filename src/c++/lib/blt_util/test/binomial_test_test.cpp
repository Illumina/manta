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

#include "blt_util/binomial_test.hh"



BOOST_AUTO_TEST_SUITE( test_binomial_test )

BOOST_AUTO_TEST_CASE( test_exact_binomial_pval )
{
    static const double tol(0.0001);

    // these tests assert a match with corresponding
    // R functions:
    const double p(0.14);
    const unsigned x(5);
    const unsigned n(12);

    BOOST_REQUIRE_CLOSE(get_binomial_twosided_exact_pval(p, x, n), 0.0361413, tol);

    BOOST_REQUIRE_CLOSE(get_binomial_gte_n_success_exact_pval(p, x, n), 0.01807065, tol);
}


BOOST_AUTO_TEST_CASE( test_simple_binomial_test )
{
    static const double alpha(0.01);
    BOOST_REQUIRE(! is_reject_binomial_twosided(alpha,0.5,1,10));
    BOOST_REQUIRE(  is_reject_binomial_twosided(alpha,0.5,10,100));

    // run the counts high enough to hit the chi-sq switchpoint:
    BOOST_REQUIRE(  is_reject_binomial_twosided(alpha,0.5,100,1000));
    BOOST_REQUIRE(! is_reject_binomial_twosided(alpha,0.1,100,1000));

    // tests to ensure that one-sided p-value from exact test are
    // working correctly.

    static const double tol(0.0001);

    {
        //simple case
        unsigned n(10);
        unsigned x(1);
        double p(0.5);

        BOOST_REQUIRE_CLOSE(get_binomial_gte_n_success_exact_pval(p, x, n),0.9990234, tol);
    }

    //simple case
    unsigned n(10);
    unsigned x(5);
    double p(0.5);

    BOOST_REQUIRE_CLOSE(get_binomial_gte_n_success_exact_pval(p, x, n),0.623046875, tol);
    BOOST_REQUIRE(! is_reject_binomial_gte_n_success_exact(0.05, p, x, n));
    BOOST_REQUIRE(  is_reject_binomial_gte_n_success_exact(0.70, p, x, n));

    // if x is 0, p-value should be 1
    BOOST_REQUIRE_CLOSE(get_binomial_gte_n_success_exact_pval(p, 0, n), 1, tol);

    // more relevant to the binomial probabilities and p-values
    // observed in somatic indel data
    n = 50;
    x = 1;
    p = 6.484e-5;

    BOOST_REQUIRE_CLOSE(get_binomial_gte_n_success_exact_pval(p, x, n), 3.23685517e-3, tol);
    BOOST_MESSAGE( "x " << x << "; n " << n << "; p" << p << " p-value is " << get_binomial_gte_n_success_exact_pval(p, x, n));
    BOOST_REQUIRE(! is_reject_binomial_gte_n_success_exact(1e-9, p, x, n));
    BOOST_REQUIRE(  is_reject_binomial_gte_n_success_exact(1e-2, p, x, n));

    x = 4;
    BOOST_REQUIRE_CLOSE(get_binomial_gte_n_success_exact_pval(p, x, n), 4.06096935e-12, tol);
    BOOST_MESSAGE( "x " << x << "; n " << n << "; p" << p << " p-value is " << get_binomial_gte_n_success_exact_pval(p, x, n));
    BOOST_REQUIRE(! is_reject_binomial_gte_n_success_exact(1e-13, p, x, n));
    BOOST_REQUIRE(  is_reject_binomial_gte_n_success_exact(1e-9,  p, x, n));
}


BOOST_AUTO_TEST_SUITE_END()
