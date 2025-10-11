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

#include <random>

template <typename T, typename RNG>
void SampleVector<T, RNG>::push(const T& val)
{
  if (_inputCount < _data.size()) {
    // initial fill of the reservoir array is deterministic:
    _data[_inputCount] = val;
  } else {
    // replace elements with gradually decreasing probability:
    std::uniform_int_distribution<unsigned> rdist(0, _inputCount);
    const unsigned                          rval(rdist(_rng));

    if (rval < _data.size()) {
      _data[rval] = val;
    }
  }
  _inputCount++;
}
