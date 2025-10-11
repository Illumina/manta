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
///

#pragma once

#include <cmath>

#include <iosfwd>
#include <limits>

/// \brief Simple on-line statistics for double values
///
/// derived From Tony Cox's IndelFinder code
///
/// TODO: there are 3 minor variants of stream_stat now... consolidate this logic via template/inheritance
///
struct stream_stat {
  // Accumulate mean and standard dev using a single pass formula
  // Uses the cancellation-friendly formulae on p.26 of
  // Higham, Accuracy & Stability of Numerical Algorithms
  // Variable names follow his
  stream_stat() : M_(0), Q_(0), max_(0), min_(0), k_(0) {}

  void reset()
  {
    M_   = 0;
    Q_   = 0;
    max_ = 0;
    min_ = 0;
    k_   = 0;
  }

  void add(const double x)
  {
    k_++;
    if (k_ == 1 || x > max_) max_ = x;
    if (k_ == 1 || x < min_) min_ = x;

    // important to do M before Q as Q uses previous iterate of M
    const double delta(x - M_);
    M_ += delta / static_cast<double>(k_);
    Q_ += delta * (x - M_);
  }

  int  size() const { return k_; }
  bool empty() const { return (k_ == 0); }

  double min() const { return ((k_ < 1) ? nan() : min_); }
  double max() const { return ((k_ < 1) ? nan() : max_); }
  double mean() const { return ((k_ < 1) ? nan() : M_); }
  double variance() const { return ((k_ < 2) ? nan() : Q_ / (static_cast<double>(k_ - 1))); }
  double sd() const { return std::sqrt(variance()); }
  double stderror() const { return sd() / std::sqrt(static_cast<double>(k_)); }

private:
  static double nan() { return std::numeric_limits<double>::quiet_NaN(); }

  double   M_;
  double   Q_;
  double   max_;
  double   min_;
  unsigned k_;
};

std::ostream& operator<<(std::ostream& os, const stream_stat& ss);
