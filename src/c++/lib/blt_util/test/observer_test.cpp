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

#include "observer.hpp"

#include <string>

BOOST_AUTO_TEST_SUITE(observer_util)

struct notifier_test : public notifier<int> {
  void increment_observers() const { notify_observers(2); }
};

struct observer_test : public observer<int> {
  observer_test() : val(0) {}

  void watch(const notifier<int>& n) { observe_notifier(n); }

private:
  void recieve_notification(const notifier<int>&, const int& n) { val += n; }

public:
  int val;
};

BOOST_AUTO_TEST_CASE(test_simple_observer)
{
  observer_test obs1;
  notifier_test not1;

  obs1.watch(not1);

  BOOST_REQUIRE_EQUAL(obs1.val, 0);
  not1.increment_observers();
  not1.increment_observers();

  BOOST_REQUIRE_EQUAL(obs1.val, 4);
}

BOOST_AUTO_TEST_CASE(test_multi_observer)
{
  observer_test obs1;
  observer_test obs2;
  notifier_test not1;
  notifier_test not2;

  obs1.watch(not1);
  obs1.watch(not2);

  obs2.watch(not1);

  BOOST_REQUIRE_EQUAL(obs1.val, 0);
  BOOST_REQUIRE_EQUAL(obs2.val, 0);

  not1.increment_observers();
  not2.increment_observers();

  BOOST_REQUIRE_EQUAL(obs1.val, 4);
  BOOST_REQUIRE_EQUAL(obs2.val, 2);
}

BOOST_AUTO_TEST_CASE(test_notifier_copy)
{
  // copy behavior for notify is designed for sane behavior in stl
  // structures, test this here
  observer_test              obs1;
  observer_test              obs2;
  std::vector<notifier_test> notvec(2);

  obs1.watch(notvec[0]);
  obs1.watch(notvec[1]);

  obs2.watch(notvec[0]);

  BOOST_REQUIRE_EQUAL(obs1.val, 0);
  BOOST_REQUIRE_EQUAL(obs2.val, 0);

  notvec[0].increment_observers();
  notvec[1].increment_observers();

  BOOST_REQUIRE_EQUAL(obs1.val, 4);
  BOOST_REQUIRE_EQUAL(obs2.val, 2);

  // simulate an alloc'd resize:
  std::vector<notifier_test> notvec2(notvec);
  notvec.clear();

  notvec2[0].increment_observers();
  notvec2[1].increment_observers();

  BOOST_REQUIRE_EQUAL(obs1.val, 8);
  BOOST_REQUIRE_EQUAL(obs2.val, 4);
}

BOOST_AUTO_TEST_SUITE_END()
