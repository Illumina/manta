// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
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

#include <cassert>

#include <algorithm>


/// return [0,1] to describe where a value is between min and max values
///
template <typename T>
struct LinearScaler
{
    LinearScaler() :
        _min(static_cast<T>(0)),
        _factor(1.)
    {}

    LinearScaler(
        const T min,
        const T max)
    {
        init(min, max);
    }

    void
    init(
        const T min,
        const T max)
    {
        assert(max>min);
        _min = min;
        _factor = (1./static_cast<double>(max-min));
    }

    double
    getScale(
        const T val) const
    {
        static const double zero(0);
        static const double one(1);
        return std::min(one, std::max(zero, static_cast<double>(val-_min)*_factor));
    }

private:
    T _min;
    double _factor;
};
