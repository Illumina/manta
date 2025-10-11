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

#include "prob_util.hpp"

BOOST_AUTO_TEST_SUITE(prob_util)

BOOST_AUTO_TEST_CASE(test_softmax)
{
  {
    static const double val(0.00001);
    const auto          tval = softMaxInverseTransform(val);
    const auto          val2 = softMaxTransform(tval);

    static const double eps = 0.00000001;
    BOOST_REQUIRE_CLOSE(val2, val, eps);
  }

  {
    static const double val(0.75);
    const auto          tval = softMaxInverseTransform(val);
    const auto          val2 = softMaxTransform(tval);

    static const double eps = 0.00000001;
    BOOST_REQUIRE_CLOSE(val2, val, eps);
  }

  {
    static const double val(500);
    static const double min(-10);
    static const double max(1000);
    const auto          tval = softMaxInverseTransform(val, min, max);
    const auto          val2 = softMaxTransform(tval, min, max);

    static const double eps = 0.00000001;
    BOOST_REQUIRE_CLOSE(val2, val, eps);
  }
}

BOOST_AUTO_TEST_CASE(test_softmax_edgecase)
{
  {
    static const double val(0.);
    const auto          tval = softMaxInverseTransform(val);
    const auto          val2 = softMaxTransform(tval);

    static const double eps = 0.00000001;
    BOOST_REQUIRE_CLOSE(val2, val, eps);
  }

  {
    static const double val(1.);
    const auto          tval = softMaxInverseTransform(val);
    const auto          val2 = softMaxTransform(tval);

    static const double eps = 0.00000001;
    BOOST_REQUIRE_CLOSE(val2, val, eps);
  }
}

BOOST_AUTO_TEST_CASE(test_softmax_factor_scale)
{
  //
  // test using softmax as a "safe" way to mult probs by interesting factors:
  //

  // positiveTest: a case where we can mult by a factor easily
  {
    static const double val(0.000001);
    static const double factor(100);
    static const double logFactor(std::log(factor));
    auto                tval = softMaxInverseTransform(val);
    tval += logFactor;
    const auto val2 = softMaxTransform(tval);

    static const double eps = 0.01;
    BOOST_REQUIRE_CLOSE(val2, (val * factor), eps);
  }

  // test2: a case where we need a "soft" max limitation to hold
  {
    static const double val(0.0049);
    static const double factor(100);
    static const double logFactor(std::log(factor));
    auto                tval = softMaxInverseTransform(val);
    tval += logFactor;
    const auto val2 = softMaxTransform(tval);

    static const double eps = 0.01;
    BOOST_REQUIRE_CLOSE(val2, 0.329944, eps);
  }
}

BOOST_AUTO_TEST_CASE(test_softmax_range)
{
  static const std::vector<double> startVal = {0.2, 0.3, 0.0001};
  auto                             val      = startVal;
  softMaxInverseRangeTransform(val.begin(), val.end());
  softMaxRangeTransform(val.begin(), val.end());

  static const double eps = 0.00000001;
  for (unsigned index(0); index < startVal.size(); ++index) {
    BOOST_REQUIRE_CLOSE(val[index], startVal[index], eps);
  }
}

BOOST_AUTO_TEST_CASE(test_prob_comp)
{
  static const float eps = 0.00000001f;

  // prob simulates a typical strong call posterior where we have strong
  // evidence for state 0 and expect a high q-score:
  const float    prob[] = {0.99999999f, 0.000000002f, 0.000000003f, 0.000000005f};
  const unsigned psize(sizeof(prob) / sizeof(float));

  const float expect = 0.00000001f;
  BOOST_REQUIRE_CLOSE(prob_comp(prob, prob + psize, 0), expect, eps);

  // uncomment this case to demo why prob_comp is used:
#if 0
    const float val1(1.-prob[0]);
    BOOST_REQUIRE_CLOSE(val1, expect, eps);
#endif
}

BOOST_AUTO_TEST_SUITE_END()
