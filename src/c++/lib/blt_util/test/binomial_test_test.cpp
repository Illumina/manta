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

#include "blt_util/binomial_test.hpp"

template <class T, size_t N>
static size_t carray_size(T (&)[N])
{
  return N;
}

BOOST_AUTO_TEST_SUITE(test_binomial_test)

BOOST_AUTO_TEST_CASE(test_exact_binomial_pval)
{
  static const double tol(0.0001);

  // these tests assert a match with corresponding
  // R functions:
  const double   p(0.14);
  const unsigned x(5);
  const unsigned n(12);

  BOOST_REQUIRE_CLOSE(get_binomial_twosided_exact_pval(p, x, n), 0.01807065, tol);

  BOOST_REQUIRE_CLOSE(get_binomial_gte_n_success_exact_pval(p, x, n), 0.01807065, tol);
}

BOOST_AUTO_TEST_CASE(test_simple_binomial_test)
{
  static const double alpha(0.01);
  BOOST_REQUIRE(!is_reject_binomial_twosided(alpha, 0.5, 1, 10));
  BOOST_REQUIRE(is_reject_binomial_twosided(alpha, 0.5, 10, 100));

  // run the counts high enough to hit the chi-sq switchpoint:
  BOOST_REQUIRE(is_reject_binomial_twosided(alpha, 0.5, 100, 1000));
  BOOST_REQUIRE(!is_reject_binomial_twosided(alpha, 0.1, 100, 1000));

  // tests to ensure that one-sided p-value from exact test are
  // working correctly.

  static const double tol(0.0001);

  {
    // simple case
    unsigned n(10);
    unsigned x(1);
    double   p(0.5);

    BOOST_REQUIRE_CLOSE(get_binomial_gte_n_success_exact_pval(p, x, n), 0.9990234, tol);
  }

  // simple case
  unsigned n(10);
  unsigned x(5);
  double   p(0.5);

  BOOST_REQUIRE_CLOSE(get_binomial_gte_n_success_exact_pval(p, x, n), 0.623046875, tol);
  BOOST_REQUIRE(!is_reject_binomial_gte_n_success_exact(0.05, p, x, n));
  BOOST_REQUIRE(is_reject_binomial_gte_n_success_exact(0.70, p, x, n));

  // if x is 0, p-value should be 1
  BOOST_REQUIRE_CLOSE(get_binomial_gte_n_success_exact_pval(p, 0, n), 1, tol);

  // more relevant to the binomial probabilities and p-values
  // observed in somatic indel data
  n = 50;
  x = 1;
  p = 6.484e-5;

  BOOST_REQUIRE_CLOSE(get_binomial_gte_n_success_exact_pval(p, x, n), 3.23685517e-3, tol);
  BOOST_TEST_MESSAGE(
      "x " << x << "; n " << n << "; p" << p << " p-value is "
           << get_binomial_gte_n_success_exact_pval(p, x, n));
  BOOST_REQUIRE(!is_reject_binomial_gte_n_success_exact(1e-9, p, x, n));
  BOOST_REQUIRE(is_reject_binomial_gte_n_success_exact(1e-2, p, x, n));

  x = 4;
  BOOST_REQUIRE_CLOSE(get_binomial_gte_n_success_exact_pval(p, x, n), 4.06096935e-12, tol);
  BOOST_TEST_MESSAGE(
      "x " << x << "; n " << n << "; p" << p << " p-value is "
           << get_binomial_gte_n_success_exact_pval(p, x, n));
  BOOST_REQUIRE(!is_reject_binomial_gte_n_success_exact(1e-13, p, x, n));
  BOOST_REQUIRE(is_reject_binomial_gte_n_success_exact(1e-9, p, x, n));
}

BOOST_AUTO_TEST_CASE(test_binomial_pvalue_0_5_manyvals)
{
  /*
   *  tests = 100
   *  trials = round(200*runif(tests))
   *  successes = round(trials*runif(tests))
   *  d = data.frame(t=seq(1, tests), trials=trials, successes=successes)
   *  d$p = 0.5
   *  d = ddply(d, .(t), function(x) {
   *      print(x)
   *      x$pval = binom.test(x$successes, x$trials, x$p, alternative="two.sided")$p.value
   *      return(x)
   *  })
   *  print(d)
   */
  //     trials successes   p         pval
  double exampledata[][4] = {
      {83, 42, 0.5, 1.000000e+00},   {104, 55, 0.5, 6.241435e-01},  {78, 38, 0.5, 9.099465e-01},
      {160, 2, 0.5, 1.762708e-44},   {168, 110, 0.5, 7.370835e-05}, {123, 24, 0.5, 5.165950e-12},
      {76, 63, 0.5, 5.044924e-09},   {95, 29, 0.5, 1.864019e-04},   {45, 29, 0.5, 7.245426e-02},
      {57, 14, 0.5, 1.538890e-04},   {48, 45, 0.5, 1.312586e-10},   {93, 58, 0.5, 2.201859e-02},
      {74, 40, 0.5, 5.613807e-01},   {193, 43, 0.5, 4.401184e-15},  {76, 13, 0.5, 5.044924e-09},
      {167, 48, 0.5, 3.825608e-08},  {171, 104, 0.5, 5.745647e-03}, {97, 42, 0.5, 2.228777e-01},
      {44, 23, 0.5, 8.803958e-01},   {2, 1, 0.5, 1.000000e+00},     {26, 25, 0.5, 8.046627e-07},
      {199, 63, 0.5, 2.466848e-07},  {152, 148, 0.5, 7.692946e-39}, {50, 42, 0.5, 1.163556e-06},
      {30, 24, 0.5, 1.430906e-03},   {173, 114, 0.5, 3.502538e-05}, {31, 1, 0.5, 2.980232e-08},
      {29, 27, 0.5, 1.624227e-06},   {20, 6, 0.5, 1.153183e-01},    {92, 5, 0.5, 2.104348e-20},
      {182, 177, 0.5, 5.284278e-46}, {14, 5, 0.5, 4.239502e-01},    {127, 46, 0.5, 2.417552e-03},
      {112, 54, 0.5, 7.769648e-01},  {146, 53, 0.5, 1.172563e-03},  {157, 40, 0.5, 5.948252e-10},
      {152, 105, 0.5, 2.899910e-06}, {196, 195, 0.5, 3.922989e-57}, {33, 12, 0.5, 1.627557e-01},
      {23, 17, 0.5, 3.468966e-02},   {104, 44, 0.5, 1.409563e-01},  {91, 87, 0.5, 2.260483e-21},
      {79, 26, 0.5, 3.183011e-03},   {16, 11, 0.5, 2.101135e-01},   {111, 63, 0.5, 1.836572e-01},
      {43, 14, 0.5, 3.153950e-02},   {157, 97, 0.5, 3.921773e-03},  {131, 60, 0.5, 3.823543e-01},
      {191, 111, 0.5, 2.968512e-02}, {99, 11, 0.5, 4.529580e-16},   {73, 58, 0.5, 4.093199e-07},
      {157, 112, 0.5, 8.922398e-08}, {51, 17, 0.5, 2.409291e-02},   {106, 29, 0.5, 3.452561e-06},
      {20, 6, 0.5, 1.153183e-01},    {187, 45, 0.5, 6.577780e-13},  {164, 149, 0.5, 6.266311e-29},
      {126, 74, 0.5, 6.093825e-02},  {163, 88, 0.5, 3.472932e-01},  {97, 46, 0.5, 6.848588e-01},
      {23, 5, 0.5, 1.062202e-02},    {126, 30, 0.5, 3.044627e-09},  {123, 42, 0.5, 5.560924e-04},
      {187, 21, 0.5, 3.631339e-29},  {174, 108, 0.5, 1.800200e-03}, {11, 2, 0.5, 6.542969e-02},
      {82, 77, 0.5, 1.204638e-17},   {11, 2, 0.5, 6.542969e-02},    {152, 138, 0.5, 8.462475e-27},
      {20, 10, 0.5, 1.000000e+00},   {139, 56, 0.5, 2.707353e-02},  {84, 59, 0.5, 2.664511e-04},
      {16, 15, 0.5, 5.187988e-04},   {87, 70, 0.5, 8.350916e-09},   {124, 55, 0.5, 2.429228e-01},
      {52, 18, 0.5, 3.648340e-02},   {136, 86, 0.5, 2.558735e-03},  {191, 59, 0.5, 1.346152e-07},
      {194, 91, 0.5, 4.297475e-01},  {84, 33, 0.5, 6.297226e-02},   {85, 84, 0.5, 4.446096e-24},
      {72, 20, 0.5, 2.077160e-04},   {191, 90, 0.5, 4.694208e-01},  {176, 169, 0.5, 2.002554e-41},
      {177, 38, 0.5, 9.936809e-15},  {112, 84, 0.5, 1.110700e-07},  {31, 26, 0.5, 1.921952e-04},
      {41, 33, 0.5, 1.122214e-04},   {182, 25, 0.5, 1.407601e-24},  {170, 3, 0.5, 1.094465e-45},
      {156, 81, 0.5, 6.890543e-01},  {11, 5, 0.5, 1.000000e+00},    {108, 12, 0.5, 1.957750e-17},
      {6, 6, 0.5, 3.125000e-02},     {116, 58, 0.5, 1.000000e+00},  {85, 40, 0.5, 6.646455e-01},
      {151, 32, 0.5, 5.591514e-13},  {170, 115, 0.5, 4.870266e-06}, {106, 10, 0.5, 8.742831e-19},
      {176, 127, 0.5, 3.617000e-09}};

  const size_t nexamples = carray_size(exampledata);
  for (size_t i = 0; i < nexamples; ++i) {
    unsigned            trials    = (unsigned)exampledata[i][0];
    unsigned            successes = (unsigned)exampledata[i][1];
    double              p         = exampledata[i][2];
    static const double tol(0.01);
    BOOST_CHECK_CLOSE(get_binomial_twosided_exact_pval(p, successes, trials), exampledata[i][3], tol);
  }
}

BOOST_AUTO_TEST_CASE(test_binomial_pvalue_many_p_manyvals)
{
  /*
   *  tests = 50
   *  trials = round(200*runif(tests))
   *  successes = round(trials*runif(tests))
   *  p = runif(tests)
   *  d = data.frame(t=seq(1, tests), trials=trials, successes=successes, p=p)
   *  d = ddply(d, .(t), function(x) {
   *      print(x)
   *      x$pval = binom.test(x$successes, x$trials, x$p, alternative="two.sided")$p.value
   *      return(x)
   *  })
   *  print(d)
   */
  //     trials successes   p         pval
  double exampledata[][4] = {
      {11, 3, 2.606046e-01, 1.000000e+00},    {119, 13, 7.146513e-01, 1.937354e-43},
      {18, 15, 9.509081e-01, 5.558266e-02},   {37, 8, 7.942178e-01, 8.049322e-14},
      {171, 140, 1.810806e-01, 3.046499e-73}, {153, 145, 2.606830e-02, 1.087423e-217},
      {105, 95, 8.408730e-02, 8.555045e-90},  {169, 125, 1.841133e-01, 1.636360e-55},
      {145, 40, 8.158598e-01, 2.012806e-45},  {167, 17, 3.623335e-01, 2.195455e-14},
      {13, 3, 9.886147e-01, 1.014544e-17},    {94, 66, 7.904077e-01, 4.206094e-02},
      {84, 70, 1.256048e-01, 4.249744e-49},   {87, 9, 7.920060e-01, 4.165875e-43},
      {151, 137, 5.341560e-01, 2.910204e-23}, {131, 5, 5.032622e-01, 6.137486e-32},
      {32, 29, 1.933186e-02, 9.393721e-47},   {154, 61, 4.758460e-02, 1.291923e-39},
      {18, 1, 8.443087e-02, 1.000000e+00},    {99, 99, 5.012314e-02, 2.012726e-129},
      {162, 4, 1.037136e-01, 2.449286e-04},   {48, 21, 8.339616e-01, 5.098986e-10},
      {195, 136, 7.164857e-01, 5.781318e-01}, {115, 98, 5.117219e-01, 2.792499e-14},
      {104, 77, 6.700700e-01, 1.442298e-01},  {68, 65, 3.900777e-01, 3.113912e-23},
      {104, 41, 6.784836e-01, 2.870279e-09},  {198, 68, 8.416943e-01, 9.580940e-56},
      {183, 17, 8.236956e-01, 1.083843e-103}, {15, 12, 8.408493e-01, 7.206539e-01},
      {81, 36, 7.994870e-01, 2.023837e-12},   {187, 40, 4.973080e-01, 1.656070e-15},
      {68, 0, 9.132837e-01, 6.177995e-73},    {65, 28, 7.287608e-01, 4.421443e-07},
      {200, 64, 6.089847e-01, 2.115927e-16},  {81, 49, 9.559001e-01, 1.785869e-22},
      {98, 21, 8.427854e-01, 4.970776e-43},   {94, 31, 4.296621e-01, 6.021201e-02},
      {118, 39, 9.823598e-01, 3.851294e-108}, {196, 121, 6.332332e-01, 6.568822e-01},
      {41, 26, 8.377987e-02, 1.804925e-18},   {7, 5, 4.714176e-01, 2.658962e-01},
      {199, 98, 2.147335e-01, 7.100855e-18},  {35, 35, 9.683290e-01, 6.288840e-01},
      {14, 6, 2.079378e-01, 5.199681e-02},    {80, 69, 3.608725e-01, 2.393712e-20},
      {71, 60, 9.829840e-02, 2.988246e-49},   {74, 48, 6.387641e-01, 9.042237e-01},
      {89, 69, 3.023493e-01, 4.819660e-20},   {7, 6, 7.700307e-01, 1.000000e+00},
  };

  const size_t nexamples = carray_size(exampledata);
  for (size_t i = 0; i < nexamples; ++i) {
    unsigned            trials    = (unsigned)exampledata[i][0];
    unsigned            successes = (unsigned)exampledata[i][1];
    double              p         = exampledata[i][2];
    static const double tol(0.01);
    BOOST_CHECK_CLOSE(get_binomial_twosided_exact_pval(p, successes, trials), exampledata[i][3], tol);
  }
}

BOOST_AUTO_TEST_CASE(test_binomial_gte_min_count)
{
  /*
  * tests <- 100

  # skew towards lower success rates, similar to what we currently see for indel error rates
  # skew the tests towards low p-values
  * test_df <- data.frame(n_trials = round(runif(tests, max = 200)),
  *                       success_rate = 10^runif(tests, min = -7, max = 0),
  *                       p_value = 10^runif(tests, min = -10, max = 0))
  * test_df$min_success <- apply(test_df, 1, function(x) qbinom(p = x["p_value"],
  *                                                             size = x["n_trials"],
  *                                                             prob = x["success_rate"],
  *                                                             lower.tail = FALSE))

  * cat(apply(test_df, 1, function(x) sprintf("{%3d, %.5e, %.5e, %3d},",
  *                                           x["n_trials"], x["success_rate"],
  *                                           x["p_value"], x["min_success"])),
  *     sep = "\n")
   */
  //  trials success_rate    p-value  min_count
  double exampledata[][4] = {{153, 2.57316e-03, 1.90497e-01, 1},  {79, 7.12531e-01, 2.08688e-04, 69},
                             {2, 2.54527e-01, 8.56125e-02, 1},    {143, 2.84603e-04, 1.13570e-07, 3},
                             {9, 7.06301e-05, 1.37532e-01, 0},    {124, 1.45277e-06, 4.07325e-10, 2},
                             {173, 1.07099e-04, 6.43000e-06, 2},  {5, 8.79078e-05, 3.78969e-07, 1},
                             {55, 7.91582e-05, 8.19433e-02, 0},   {99, 7.32748e-04, 8.09058e-02, 0},
                             {7, 9.42784e-02, 2.01747e-04, 4},    {129, 1.50969e-04, 1.04377e-10, 4},
                             {121, 5.45886e-01, 7.98495e-09, 96}, {188, 1.43810e-03, 4.87911e-01, 0},
                             {114, 1.78681e-01, 1.79785e-03, 33}, {130, 1.39812e-07, 3.98930e-03, 0},
                             {97, 5.02856e-03, 3.02772e-07, 7},   {183, 1.43512e-05, 8.16882e-09, 2},
                             {185, 2.90423e-03, 2.18749e-05, 5},  {84, 2.50545e-03, 5.12018e-03, 2},
                             {116, 9.91294e-04, 2.63842e-05, 3},  {109, 8.89566e-06, 3.43128e-07, 2},
                             {39, 7.61711e-05, 3.75268e-07, 2},   {193, 3.68159e-03, 1.62027e-08, 9},
                             {142, 2.59871e-06, 1.50306e-05, 1},  {182, 4.02888e-04, 9.75801e-01, 0},
                             {141, 2.14857e-07, 6.16874e-06, 1},  {15, 9.18554e-07, 4.03038e-05, 0},
                             {149, 2.28428e-07, 2.12637e-04, 0},  {191, 2.76298e-07, 6.71069e-09, 1},
                             {107, 1.31901e-04, 1.04454e-09, 4},  {65, 1.85190e-07, 5.78113e-02, 0},
                             {141, 2.66733e-06, 1.44922e-02, 0},  {190, 1.03013e-05, 1.53499e-09, 2},
                             {49, 3.26117e-05, 1.84198e-09, 2},   {127, 1.25834e-02, 2.48552e-05, 9},
                             {42, 4.28448e-07, 7.31442e-08, 1},   {60, 6.90027e-03, 1.05096e-07, 7},
                             {144, 2.90896e-02, 6.76604e-06, 15}, {138, 4.18742e-01, 4.95768e-02, 67},
                             {117, 9.73191e-06, 4.88660e-04, 1},  {24, 6.70017e-07, 2.68855e-09, 1},
                             {102, 1.58112e-05, 8.01416e-02, 0},  {192, 1.25658e-02, 2.44385e-09, 16},
                             {50, 3.17352e-03, 7.21933e-09, 6},   {27, 1.00049e-04, 8.91732e-06, 1},
                             {47, 6.11921e-04, 6.35839e-03, 1},   {168, 2.29714e-03, 6.40052e-07, 6},
                             {73, 1.26782e-06, 3.10391e-10, 2},   {195, 1.88220e-04, 8.98213e-03, 1},
                             {1, 1.92718e-07, 3.95304e-06, 0},    {155, 7.96789e-06, 5.51888e-01, 0},
                             {20, 1.71476e-01, 5.86904e-10, 16},  {1, 2.43012e-04, 1.06529e-05, 1},
                             {149, 1.07803e-07, 7.93280e-09, 1},  {11, 4.73269e-05, 1.00482e-10, 2},
                             {72, 8.34266e-01, 6.74543e-09, 72},  {198, 4.17837e-05, 8.69572e-02, 0},
                             {142, 3.85465e-07, 6.11033e-04, 0},  {105, 1.42848e-02, 1.19440e-08, 12},
                             {161, 3.83509e-04, 7.23564e-03, 1},  {103, 9.24136e-04, 3.22721e-08, 5},
                             {107, 1.13859e-05, 9.39218e-01, 0},  {4, 5.05847e-07, 1.06546e-01, 0},
                             {184, 6.77215e-03, 2.78220e-10, 13}, {142, 1.20315e-01, 7.06077e-06, 36},
                             {113, 1.53915e-02, 6.81864e-07, 11}, {151, 5.17850e-07, 8.65573e-04, 0},
                             {87, 3.40905e-05, 6.05011e-09, 2},   {106, 8.37584e-03, 2.17082e-02, 3},
                             {14, 7.23873e-07, 2.66869e-04, 0},   {40, 1.44629e-06, 1.74549e-04, 0},
                             {196, 2.15504e-06, 4.37458e-06, 1},  {75, 1.74346e-02, 1.41908e-04, 7},
                             {14, 5.88490e-04, 1.12071e-04, 1},   {175, 9.32951e-03, 1.76643e-07, 11},
                             {39, 1.99337e-03, 3.16462e-05, 3},   {15, 5.46313e-01, 1.00483e-04, 15},
                             {24, 6.67693e-01, 3.52558e-04, 23},  {129, 2.56853e-03, 4.44043e-04, 3},
                             {99, 4.41006e-05, 2.10157e-01, 0},   {185, 4.27154e-04, 8.46782e-03, 1},
                             {39, 6.18477e-03, 2.52357e-02, 1},   {167, 3.26215e-02, 1.57803e-10, 25},
                             {101, 5.79464e-04, 5.95059e-09, 4},  {13, 9.57207e-07, 1.10881e-03, 0},
                             {124, 1.96723e-04, 7.19052e-02, 0},  {149, 3.75192e-04, 2.53588e-06, 3},
                             {103, 8.99213e-06, 3.28940e-03, 0},  {121, 4.10346e-03, 9.68812e-01, 0},
                             {196, 7.97230e-07, 3.46923e-01, 0},  {9, 5.55411e-04, 9.63213e-07, 2},
                             {134, 2.16316e-04, 2.82896e-07, 3},  {99, 1.03844e-05, 2.72108e-06, 1},
                             {140, 1.71302e-03, 2.48607e-01, 0},  {196, 1.60494e-04, 7.52108e-10, 4},
                             {147, 1.10965e-05, 4.24813e-01, 0},  {34, 5.76217e-03, 7.92384e-09, 6},
                             {153, 1.05309e-07, 2.31016e-01, 0},  {112, 1.39439e-04, 5.75978e-09, 3}};

  const size_t nexamples = carray_size(exampledata);
  for (size_t i = 0; i < nexamples; ++i) {
    unsigned trials       = (unsigned)exampledata[i][0];
    double   success_rate = exampledata[i][1];
    double   p_val        = exampledata[i][2];
    BOOST_CHECK_EQUAL(
        (unsigned)min_count_binomial_gte_exact(p_val, success_rate, trials), 1 + exampledata[i][3]);
  }
}
BOOST_AUTO_TEST_SUITE_END()
