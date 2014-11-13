// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
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
struct SampleVector
{
    /// \param initrng c++11 <random> rng generator, see std::shuffle for detailed doc of similar parameter
    SampleVector(
        const unsigned initSize,
        RNG& initRng)
        : _inputCount(0),
          _data(initSize,0),
          _rng(initRng)
    {}

    void
    push(const T& val);

    const std::vector<T>&
    data() const
    {
        return _data;
    }

private:
    unsigned _inputCount;
    std::vector<T> _data;
    RNG& _rng;
};

#include "SampleVectorImpl.hh"
