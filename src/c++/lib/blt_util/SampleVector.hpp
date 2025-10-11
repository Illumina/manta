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

#include <vector>

/// random sub-sampling array
///
/// This is an array with sized fixed to S at instantiation time. The
/// array accepts N (N>=S) objects as input. For any N>=S, the array
/// contains a given input object with probability S/N.
///
/// This behavior is implemented via standard reservoir sampling
///
template <typename T, typename RNG>
struct SampleVector {
  /// \param initrng c++11 <random> rng generator, see std::shuffle for detailed doc of similar parameter
  SampleVector(const unsigned initSize, RNG& initRng) : _inputCount(0), _data(initSize, 0), _rng(initRng) {}

  void push(const T& val);

  const std::vector<T>& data() const { return _data; }

private:
  unsigned       _inputCount;
  std::vector<T> _data;
  RNG&           _rng;
};

#include "SampleVectorImpl.hpp"
