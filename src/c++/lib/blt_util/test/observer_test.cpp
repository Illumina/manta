//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
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
