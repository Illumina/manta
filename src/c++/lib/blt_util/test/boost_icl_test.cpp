//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
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

#include "boost/test/unit_test.hpp"

#include "boost/icl/discrete_interval.hpp"
#include "boost/icl/interval_map.hpp"

#include <set>
#include <string>

using namespace boost::icl;

BOOST_AUTO_TEST_SUITE(test_boost_icl)

BOOST_AUTO_TEST_CASE(boost_icl_test_intersect)
{
  discrete_interval<int> interval1 = interval<int>::right_open(3, 7);
  discrete_interval<int> interval2 = interval<int>::right_open(7, 8);
  discrete_interval<int> interval3 = interval<int>::right_open(6, 8);

  BOOST_REQUIRE(!intersects(interval1, interval2));
  BOOST_REQUIRE(intersects(interval1, interval3));
}

BOOST_AUTO_TEST_CASE(boost_icl_test_map)
{
  typedef std::set<std::string>       test_val;
  typedef interval_map<int, test_val> map_t;
  map_t                               test_map;
  test_val                            foo, bar;
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
            for (const std::string& s : begin->second)
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
            for (const std::string& s : begin->second)
            {
                std::cerr << " " << s;
            }
            std::cerr << "\n";
        }
    }
#endif
}

BOOST_AUTO_TEST_CASE(boost_icl_test_map2)
{
  typedef interval_map<int, unsigned> map_t;
  map_t                               test_map;

  test_map.add(std::make_pair(interval<int>::right_open(3, 7), 1u));
  test_map.add(std::make_pair(interval<int>::right_open(4, 5), 2u));
  test_map.add(std::make_pair(interval<int>::right_open(6, 8), 1u));

  BOOST_REQUIRE_EQUAL(test_map.iterative_size(), 5u);

  auto b(test_map.begin());
  BOOST_REQUIRE_EQUAL(b->first.lower(), 3);
  BOOST_REQUIRE_EQUAL(b->first.upper(), 4);
  BOOST_REQUIRE_EQUAL(b->second, 1u);
  ++b;
  BOOST_REQUIRE_EQUAL(b->first.lower(), 4);
  BOOST_REQUIRE_EQUAL(b->first.upper(), 5);
  BOOST_REQUIRE_EQUAL(b->second, 3u);
  ++b;
  BOOST_REQUIRE_EQUAL(b->first.lower(), 5);
  BOOST_REQUIRE_EQUAL(b->first.upper(), 6);
  BOOST_REQUIRE_EQUAL(b->second, 1u);
  ++b;
  BOOST_REQUIRE_EQUAL(b->first.lower(), 6);
  BOOST_REQUIRE_EQUAL(b->first.upper(), 7);
  BOOST_REQUIRE_EQUAL(b->second, 2u);
  ++b;
  BOOST_REQUIRE_EQUAL(b->first.lower(), 7);
  BOOST_REQUIRE_EQUAL(b->first.upper(), 8);
  BOOST_REQUIRE_EQUAL(b->second, 1u);
}

BOOST_AUTO_TEST_CASE(boost_icl_test_map3)
{
  // recreates ICL bug reported for OS X 10.9/clang 3.5 and boost 1.53
  // (probably) related: https://github.com/Astron/Astron/issues/213
  // intent is to test for an assertion from boost ICL, boost value
  // tests below are just placeholders.
  //
  // known to pass with boost 1.56+ and fail in boost 1.55-
  //

  typedef interval_map<int, unsigned> map_t;
  map_t                               test_map;

  test_map.add(std::make_pair(interval<int>::right_open(3, 5), 1u));
  test_map.add(std::make_pair(interval<int>::right_open(8, 9), 2u));
  test_map.add(std::make_pair(interval<int>::right_open(1, 12), 1u));

  BOOST_REQUIRE_EQUAL(test_map.iterative_size(), 5u);
}

BOOST_AUTO_TEST_SUITE_END()
