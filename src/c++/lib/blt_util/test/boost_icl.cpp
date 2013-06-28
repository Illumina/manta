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

#include "boost/foreach.hpp"
#include "boost/icl/discrete_interval.hpp"
#include "boost/icl/interval_map.hpp"

#include <set>
#include <string>

using namespace boost::icl;

//
// this is more of a test/demo for boost:icl then an actual test....
//

BOOST_AUTO_TEST_SUITE( test_boost_icl )


BOOST_AUTO_TEST_CASE( boost_icl_test_intersect )
{
    discrete_interval<int> interval1 = interval<int>::right_open(3,7);
    discrete_interval<int> interval2 = interval<int>::right_open(7,8);
    discrete_interval<int> interval3 = interval<int>::right_open(6,8);

    BOOST_REQUIRE(! intersects(interval1,interval2));
    BOOST_REQUIRE( intersects(interval1,interval3));
}


BOOST_AUTO_TEST_CASE( boost_icl_test_map )
{
    typedef std::set<std::string> test_val;
    typedef interval_map<int,test_val> map_t;
    map_t test_map;
    test_val foo,bar;
    foo.insert("foo");
    bar.insert("bar");


    // commenting out this code so that we can compile with clang, apparently
    // clang and boost::icl disagree -- not worth fixing this uniless icl
    // is used in a production code path
    //

#if 0
    test_map.add(std::make_pair(interval<int>::right_open(3,7), foo));
    test_map.add(std::make_pair(interval<int>::right_open(6,8), bar));

    {
        map_t::const_iterator begin(test_map.find(5));
        map_t::const_iterator end(test_map.find(7));
        for (; begin!=end; ++begin)
        {
            std::cerr << "CCC" << begin->first;
            BOOST_FOREACH(const std::string& s, begin->second)
            {
                std::cerr << " " << s;
            }
            std::cerr << "\n";
        }
    }

    test_map.erase(std::make_pair(interval<int>::right_open(6,8), bar));
    {
        map_t::const_iterator begin(test_map.find(5));
        map_t::const_iterator end(test_map.find(7));
        for (; begin!=end; ++begin)
        {
            std::cerr << "CCC" << begin->first;
            BOOST_FOREACH(const std::string& s, begin->second)
            {
                std::cerr << " " << s;
            }
            std::cerr << "\n";
        }
    }
#endif
    }


BOOST_AUTO_TEST_SUITE_END()

