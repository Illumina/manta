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

/// \file
/// \author Chris Saunders
/// \author Mitch Bekritsky
///

#pragma once

/// \brief Two-sided binomial exact probability
///
/// This is a two sided binomial exact pval wherein we find the
/// prob of n_success or more extreme number of successes and then
/// double it.
///
///
double get_binomial_twosided_exact_pval(const double p, const unsigned n_success, const unsigned n_trials);

/// \brief Two-sided binomial exact test
///
bool is_reject_binomial_twosided_exact(
    const double alpha, const double p, const unsigned n_success, const unsigned n_trials);

/// \brief Two-sided binomial test chi-sqr approximation
bool is_reject_binomial_twosided_chi_sqr(
    const double alpha, const double p, const unsigned n_success, const unsigned n_trials);

/// \brief Two-sided binomial test
///
/// Find the probability of n_success or more extreme success under B(n_trial,p)
///
/// This function chooses from the two testing methods above (exact/approx) based on trial size
///
bool is_reject_binomial_twosided(
    const double alpha, const double p, const unsigned n_success, const unsigned n_trials);

/// \brief One-sided binomial exact probability
///
/// Return probability of \p n_success or more given B(n_trials,p)
///
/// This should match R code:
/// \code{.R}
/// pbinom((n_success-1),n_trials,p,lower.tail=FALSE)
/// \endcode
double get_binomial_gte_n_success_exact_pval(
    const double p, const unsigned n_success, const unsigned n_trials);

/// \brief One-sided binomial exact test
///
/// Test whether \p n_success or greater can be rejected under
/// a null hypothesis of B(n_trials,p)
///
/// This should match R code:
/// \code{.R}
/// binom.test(n_success, n_trials, p, "greater")$p.value <= alpha
/// \endcode
bool is_reject_binomial_gte_n_success_exact(
    const double alpha, const double p, const unsigned n_success, const unsigned n_trials);

/// Return the minimum number of successes to reject the null hypothesis
/// with a p-value of at most \p alpha for a given error rate and number of trials
///
/// This should match R code:
/// \code{.R}
/// 1 + qbinom(alpha, n_trials, p, lower.tail = FALSE)
/// \endcode
double min_count_binomial_gte_exact(const double alpha, const double p, const unsigned n_trials);
