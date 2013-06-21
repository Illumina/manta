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

#include "observer.hh"

#include <string>

BOOST_AUTO_TEST_SUITE( observer_util )


struct notifier_test : public notifier<int>
{
    void
    increment_observers() const
    {
        notify_observers(2);
    }
};



struct observer_test : public observer<int>
{
    observer_test() :
        val(0)
    {}

    void
    watch(const notifier<int>& n)
    {
        observe_notifier(n);
    }

private:
    void
    recieve_notification(const notifier<int>&,
                         const int& n) {
        val += n;
    }

public:
    int val;
};




BOOST_AUTO_TEST_CASE( test_simple_observer )
{
    observer_test obs1;
    notifier_test not1;

    obs1.watch(not1);

    BOOST_REQUIRE_EQUAL(obs1.val,0);
    not1.increment_observers();
    not1.increment_observers();

    BOOST_REQUIRE_EQUAL(obs1.val,4);
}

#if 0
BOOST_AUTO_TEST_CASE( test_multi_observer )
{
    observer_test obs1;
    observer_test obs2;
    notifier_test not1;
    notifier_test not2;

    obs1.watch(not1);
    obs1.watch(not2);

    obs2.watch(not1);

    BOOST_REQUIRE_EQUAL(obs1.val,0);
    BOOST_REQUIRE_EQUAL(obs2.val,0);

    not1.increment_observers();
    not2.increment_observers();

    BOOST_REQUIRE_EQUAL(obs1.val,4);
    BOOST_REQUIRE_EQUAL(obs2.val,2);
}
#endif

BOOST_AUTO_TEST_SUITE_END()

