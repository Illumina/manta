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

#include "stage_manager.hpp"

//#define DEBUG_SM_TEST

#ifdef DEBUG_SM_TEST
#include <iostream>

namespace {
std::ostream& log_os(std::cerr);
}
#endif

BOOST_AUTO_TEST_SUITE(test_stage_manager)

// create a standard stage size arrangement for testing:
//
// This returns a tree with:
//   stage 1 following 10 bases behind the root
//   stage 2 following 20 bases behind stage 1
//   stage 3 following 20 bases behind the root
//
static stage_data get_test_stage_data()
{
  stage_data sd;
  sd.add_stage(0);
  sd.add_stage(1, 0, 10);
  sd.add_stage(2, 1, 20);
  sd.add_stage(3, 0, 20);

  return sd;
}

static void stage_data_shift_test(const stage_data& sd, const int stage_id, const unsigned expect)
{
  const unsigned result(sd.get_stage_id_shift(stage_id));
  BOOST_CHECK_EQUAL(result, expect);
}

BOOST_AUTO_TEST_CASE(test_stage_data_dist)
{
  const stage_data sd(get_test_stage_data());

  stage_data_shift_test(sd, 0, 0);
  stage_data_shift_test(sd, 1, 10);
  stage_data_shift_test(sd, 2, 30);
  stage_data_shift_test(sd, 3, 20);
}

BOOST_AUTO_TEST_CASE(test_stage_data_bad_parent)
{
  stage_data sd;
  BOOST_CHECK_THROW(sd.add_stage(1, 0, 10), std::exception);
}

BOOST_AUTO_TEST_CASE(test_stage_data_bad_id)
{
  stage_data sd;
  sd.add_stage(1);

  BOOST_CHECK_THROW(sd.get_stage_id_shift(0), std::exception);
  BOOST_CHECK_NO_THROW(sd.get_stage_id_shift(1));
  BOOST_CHECK_THROW(sd.get_stage_id_shift(2), std::exception);
}

/// \brief Minimal pos_processor object used to test stage manager
///
/// Note that this object is itself part of the test infrastructure by
/// asserting:
///   1. ...that pos increases for each stage
///   2. TODO: ..the expected relationship (stage-to-root distance vs expect)
///              of all stages as process_pos is called
///
struct test_pos_processor : public pos_processor_base {
  //
  // TODO: finish setting up stage relationship checking...
  //
  // pos_processor wouldn't normally need this info, but we use
  // it to test expected stage position relationships
  //
  // test_pos_processor(const stage_data& sd, const pos_range& pr)
  //

  void process_pos(const int stage_no, const pos_t pos)
  {
#ifdef DEBUG_SM_TEST
    log_os << "process_pos stage_no: " << stage_no << " pos: " << pos << "\n";
#endif

    // assert that pos for each stage does not repeat or decrease:
    spos_t::const_iterator i(stage_pos.find(stage_no));
    if (i != stage_pos.end()) {
      BOOST_CHECK(pos > (i->second));
    }
    stage_pos[stage_no] = pos;
  }

  typedef std::map<int, pos_t> spos_t;
  spos_t                       stage_pos;
};

BOOST_AUTO_TEST_CASE(test_stage_manager)
{
  const stage_data   sd(get_test_stage_data());
  const pos_range    report_range(0, 60);
  test_pos_processor tpp;

  stage_manager sman(sd, report_range, tpp);

  sman.handle_new_pos_value(40);

  BOOST_CHECK_EQUAL(tpp.stage_pos[0], 40);
  BOOST_CHECK_EQUAL(tpp.stage_pos[1], 30);
  BOOST_CHECK_EQUAL(tpp.stage_pos[2], 10);
  BOOST_CHECK_EQUAL(tpp.stage_pos[3], 20);
}

BOOST_AUTO_TEST_CASE(test_stage_manager_reset)
{
  const stage_data   sd(get_test_stage_data());
  const pos_range    report_range(0, 60);
  test_pos_processor tpp;

  stage_manager sman(sd, report_range, tpp);

  sman.reset();

  for (int i(0); i < 4; ++i) {
    BOOST_CHECK_EQUAL(tpp.stage_pos[i], 59);
  }
}

BOOST_AUTO_TEST_SUITE_END()
