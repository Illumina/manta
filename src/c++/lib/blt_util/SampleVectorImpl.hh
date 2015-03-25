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

#include <random>



template <typename T, typename RNG>
void
SampleVector<T,RNG>::
push(const T& val)
{
    if (_inputCount < _data.size())
    {
        // initial fill of the reservoir array is deterministic:
        _data[_inputCount] = val;
    }
    else
    {
        // replace elements with gradually decreasing probability:
        std::uniform_int_distribution<unsigned> rdist(0,_inputCount);
        const unsigned rval(rdist(_rng));

        if (rval < _data.size())
        {
            _data[rval] = val;
        }
    }
    _inputCount++;
}
